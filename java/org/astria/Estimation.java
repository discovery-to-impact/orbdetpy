/*
 * Estimation.java - Implementation of estimation algorithms.
 * Copyright (C) 2018-2019 University of Texas
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.astria;

import com.google.gson.GsonBuilder;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import org.astria.DataManager;
import org.astria.ManualPropagation;
import org.astria.Measurements;
import org.astria.PropagatorBuilder;
import org.astria.Settings;
import org.astria.CAR;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.DiagonalMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.LUDecomposition;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.sequential.ConstantProcessNoise;
import org.orekit.estimation.sequential.KalmanEstimation;
import org.orekit.estimation.sequential.KalmanEstimator;
import org.orekit.estimation.sequential.KalmanEstimatorBuilder;
import org.orekit.estimation.sequential.KalmanObserver;
import org.hipparchus.distribution.continuous.ChiSquaredDistribution;
import org.orekit.forces.ForceModel;
import org.orekit.forces.gravity.NewtonianAttraction;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.conversion.DormandPrince853IntegratorBuilder;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.Constants;
import org.orekit.utils.ParameterDriver;
import org.orekit.utils.ParameterDriversList;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;
import org.orekit.frames.Frame;
import org.orekit.errors.OrekitException;

import java.io.FileWriter;
import java.io.IOException;


public class Estimation
{
    protected Settings odcfg; //eventually delete
    protected Measurements odobs;

    protected String[] meanames;
    protected boolean combmeas;

    protected JSONResults results;
    protected ArrayList<JSONResults> resultsArr = new ArrayList<JSONResults>();;
	protected ArrayList<Measurements.JSONMeasurement> unassociatedObsJSON;
	protected ArrayList<Measurements.JSONMeasurement> rawMeasurements;
    

    public final static String DMC_ACC_ESTM = "DMCaccest";
    public final static String DMC_ACC_PROP = "DMCaccprop";

    
	ArrayList<ArrayList<SCstate>> states = new ArrayList<ArrayList<SCstate>>();
	ArrayList<SCstate> promotedTracks = new ArrayList<SCstate>();
    
    
    public Estimation(String[] cfgjson, String obsjson)
    {

    for(int i = 0; i < cfgjson.length; i++)
    {

		states.add(i, new ArrayList<SCstate>());

		states.get(i).add(new SCstate());
		
		states.get(i).get(0).odcfg = Settings.loadJSON(cfgjson[i]);

		if (states.get(i).get(0).odcfg.Estimation.Filter == null)
			states.get(i).get(0).odcfg.Estimation.Filter = "UKF";

		if (states.get(i).get(0).odcfg.Estimation.Filter.equals("UKF") && states.get(i).get(0).odcfg.Gravity.Degree >= 2 && states.get(i).get(0).odcfg.Gravity.Order >= 0)
			states.get(i).get(0).odcfg.forces.add(0, new NewtonianAttraction(Constants.EGM96_EARTH_MU));

		if (states.get(i).get(0).odcfg.Estimation.NoiseTimeDelta <= 0.0)
			states.get(i).get(0).odcfg.Estimation.NoiseTimeDelta = 10.0;

    }
		odobs = Measurements.loadJSON(states.get(0).get(0).odcfg, obsjson);

		
		// Remove because these measurements are not considered when estimating
		// If measurements are used it is for CAR, which only needs the level of noise in the measurements.
		// This assumes equal noise for all stations.
		Map<String, Settings.JSONMeasurement> AllMeasurements = states.get(0).get(0).odcfg.Measurements;
		
		for(int i = 0; i < states.size(); i++)
		{
			AllMeasurements = states.get(i).get(0).odcfg.Measurements;

			
			if(AllMeasurements.get("RightAscensionRate")!= null)
			{
				states.get(i).get(0).RaNoise = AllMeasurements.remove("RightAscensionRate").Error[0]; 
			}
			
			if(AllMeasurements.get("DeclinationRate")!= null)
			{
				states.get(i).get(0).DecNoise = AllMeasurements.remove("DeclinationRate").Error[0]; 
			}
		}

		    
		
		meanames = AllMeasurements.keySet().toArray(new String[0]);
				
		combmeas = meanames[0].equals("Azimuth") || meanames[0].equals("Elevation") ||
		    meanames[0].equals("RightAscension") || meanames[0].equals("Declination") ||
		    meanames[0].equals("PositionVelocity");

    }

    public String[] determineOrbit()
    {

	results = new JSONResults();
	results.Filter = states.get(0).get(0).odcfg.Estimation.Filter;


	
	if (results.Filter.equals("UKF"))
	    new UnscentedKalmanFilter().determineOrbit();
	else
	    new ExtendedKalmanFilter().determineOrbit();


	String[] JSONOutput = new String[2*resultsArr.size()+1];

	for(int i = 0; i < resultsArr.size(); i++)
	{
		JSONOutput[i] =new GsonBuilder().setPrettyPrinting().create().toJson(resultsArr.get(i));
		JSONOutput[i+resultsArr.size()] =new GsonBuilder().setPrettyPrinting().create().toJson(promotedTracks.get(i).associatedObsJSON);

	}
	
	JSONOutput[2*resultsArr.size()] = new GsonBuilder().setPrettyPrinting().create().toJson(rawMeasurements);
		
	return(JSONOutput);
    }

    
    
    
    
    
    protected class ExtendedKalmanFilter implements KalmanObserver
    {
	protected void determineOrbit()
	{
	    double[] Xi = odcfg.getInitialState();
	    CartesianOrbit X0 = new CartesianOrbit(new PVCoordinates(
						       new Vector3D(Xi[0], Xi[1], Xi[2]),
						       new Vector3D(Xi[3], Xi[4], Xi[5])),
				   		odcfg.propframe, new AbsoluteDate(
						       DateTimeComponents.parseDateTime(odcfg.Propagation.Start),
						       DataManager.utcscale), Constants.EGM96_EARTH_MU);

	    PropagatorBuilder prop = new PropagatorBuilder(odcfg, X0, new DormandPrince853IntegratorBuilder(
							       odcfg.Integration.MinTimeStep, odcfg.Integration.MaxTimeStep, 1.0),
							   PositionAngle.MEAN, 10.0);
	    prop.setMass(odcfg.SpaceObject.Mass);
	    for (ForceModel fm : odcfg.forces)
		prop.addForceModel(fm);

	    ParameterDriversList plst = prop.getPropagationParametersDrivers();
	    for (Settings.EstimatedParameter ep : odcfg.estparams)
	    {
		ParameterDriver pdrv = new ParameterDriver(ep.name, ep.value, 1.0, ep.min, ep.max);
		pdrv.setSelected(true);
		plst.add(pdrv);
	    }

	    AttitudeProvider attprov = odcfg.getAttitudeProvider();
	    if (attprov != null)
		prop.setAttitudeProvider(attprov);

	    KalmanEstimatorBuilder build = new KalmanEstimatorBuilder();
	    build.addPropagationConfiguration(prop, new ConstantProcessNoise(new DiagonalMatrix(odcfg.Estimation.Covariance),
									     odcfg.getProcessNoiseMatrix()));

	    KalmanEstimator filt = build.build();
	    filt.setObserver(this);
	    NumericalPropagator est = filt.processMeasurements(odobs.measobjs)[0];

	    results.Propagation.Time = odcfg.Propagation.End;
	    results.Propagation.State = new double[odcfg.estparams.size() + 6];

	    SpacecraftState ssta = est.propagate(new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.End),
								  DataManager.utcscale));
	    if (ssta.getA() <= Constants.WGS84_EARTH_EQUATORIAL_RADIUS)
		throw(new RuntimeException(String.format("Invalid semi-major axis %f", ssta.getA())));

	    PVCoordinates pvc = ssta.getPVCoordinates(odcfg.propframe);
	    System.arraycopy(pvc.getPosition().toArray(), 0, results.Propagation.State, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, results.Propagation.State, 3, 3);
	    if (odcfg.estparams.size() > 0)
		System.arraycopy(results.Estimation.get(results.Estimation.size() - 1).EstimatedState, 6,
		results.Propagation.State, 6, odcfg.estparams.size());
	}

	public void evaluationPerformed(KalmanEstimation est)
	{
	    int n = est.getCurrentMeasurementNumber() - 1;
	    if (!combmeas)
		n /= meanames.length;

	    String k;
	    JSONResults.JSONEstimation res;
	    if (results.Estimation.size() <= n)
	    {
		k = meanames[0];
		res = results.new JSONEstimation();
		res.Time = odobs.rawmeas[n].Time;
		results.Estimation.add(res);
	    }
	    else
	    {
		k = meanames[1];
		res = results.Estimation.get(n);
	    }

	    SpacecraftState ssta = est.getPredictedSpacecraftStates()[0];
	    if (ssta.getA() <= Constants.WGS84_EARTH_EQUATORIAL_RADIUS)
		throw(new RuntimeException(String.format("Invalid semi-major axis %f", ssta.getA())));

	    PVCoordinates pvc = ssta.getPVCoordinates();
	    res.EstimatedState = new double[odcfg.estparams.size() + 6];
	    System.arraycopy(pvc.getPosition().toArray(), 0, res.EstimatedState, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, res.EstimatedState, 3, 3);

	    res.EstimatedAcceleration = new double[3];
	    System.arraycopy(pvc.getAcceleration().toArray(), 0, res.EstimatedAcceleration, 0, 3);

	    int i = 6;
	    List<ParameterDriversList.DelegatingDriver> plst = est.getEstimatedPropagationParameters().getDrivers();
	    for (Settings.EstimatedParameter ep : odcfg.estparams)
		for (ParameterDriversList.DelegatingDriver dd : plst)
		    if (dd.getName().equals(ep.name))
			res.EstimatedState[i++] = dd.getValue();

	    if (ssta.hasAdditionalState(Estimation.DMC_ACC_PROP))
	    {
		double[] accric = ssta.getAdditionalState(Estimation.DMC_ACC_PROP);
		System.arraycopy(accric, 0, res.EstimatedState, odcfg.estparams.size() + 3, 3);
	    }

	    double[] pre = est.getPredictedMeasurement().getEstimatedValue();
	    double[] pos = est.getCorrectedMeasurement().getEstimatedValue();

	    if (combmeas)
	    {
		for (i = 0; i < meanames.length; i++)
		{
		    if (meanames.length == 1)
		    {
			res.PreFit.put(meanames[i], pre);
			res.PostFit.put(meanames[i], pos);
		    }
		    else
		    {
			res.PreFit.put(meanames[i], new double[] {pre[i]});
			res.PostFit.put(meanames[i], new double[] {pos[i]});
		    }
		}

		res.InnovationCovariance = est.getPhysicalInnovationCovarianceMatrix().getData();
	    }
	    else
	    {
		res.PreFit.put(k, pre);
		res.PostFit.put(k, pos);

		if (res.InnovationCovariance == null)
		    res.InnovationCovariance = new double[][]{{
			    est.getPhysicalInnovationCovarianceMatrix().getData()[0][0], 0.0}, {0.0, 0.0}};
		else
		    res.InnovationCovariance[1][1] = est.getPhysicalInnovationCovarianceMatrix().getData()[0][0];
	    }

	    res.EstimatedCovariance = est.getPhysicalEstimatedCovarianceMatrix().getData();
	}
    }

    protected class UnscentedKalmanFilter
    {
	protected void determineOrbit()
	{
		try {
			FileWriter writer = new FileWriter("../../../Matlab/outputStates.txt"); 
			writer.close();
			writer = new FileWriter("../../../Matlab/MarginalEventPsi.txt"); 	
			writer.close();
			writer = new FileWriter("../../../Matlab/SumJPDALiklihood.txt"); 
			writer.close();
			writer = new FileWriter("../../../Matlab/betaSatMeas.txt"); 
			writer.close();
			writer = new FileWriter("../../../Matlab/HypWeight.txt"); 
			writer.close();



		} catch (IOException e) {
			
		e.printStackTrace();
	
		} 
		
		
		int numsta = states.get(0).get(0).odcfg.estparams.size() + 6;
	    int numsig = 2*numsta;
	    int veclen = numsta*numsig;
		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			SCstate singleHyp = states.get(objNum).get(0);

			
			Settings odcfg = singleHyp.odcfg;
			singleHyp.P = new DiagonalMatrix(odcfg.Estimation.Covariance);
			singleHyp.Q = odcfg.getProcessNoiseMatrix();
			
			singleHyp.Rsize = 0;
		    for (String s: meanames)
		    	singleHyp.Rsize += odcfg.Measurements.get(s).Error.length;

		    singleHyp.R = new Array2DRowRealMatrix(singleHyp.Rsize, singleHyp.Rsize);
		    for (int i = 0, j = 0; i < meanames.length; i++)
		    {
			Settings.JSONMeasurement jm = odcfg.Measurements.get(meanames[i]);
			for (int k = 0; k < jm.Error.length; k++)
			{
				singleHyp.R.setEntry(j, j, jm.Error[k]*jm.Error[k]);
			    j++;
			}
		    }

		    singleHyp.Xi = odcfg.getInitialState();
		    singleHyp.epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.Start),
							  DataManager.utcscale);
	
		    //if (odobs.rawmeas.length > 0)
		    //{
			//AbsoluteDate tmto = new AbsoluteDate(DateTimeComponents.parseDateTime(odobs.rawmeas[0].Time),
			//				     DataManager.utcscale);
			//if (Math.abs(tmto.durationFrom(states.get(objNum).epoch)) > 1.0)
			//{
			    //ManualPropagation prop0 = new ManualPropagation(odcfg, 6);
			    //double[] Xout = prop0.propagate(0, Arrays.copyOfRange(states.get(objNum).Xi, 0, 6), tmto.durationFrom(states.get(objNum).epoch)); //what about propagating covar also?
			    //states.get(objNum).epoch = new AbsoluteDate(tmto, 0.0);
			    //System.arraycopy(Xout, 0, states.get(objNum).Xi, 0, 6);
			//}
		    //}
	


		    singleHyp.xhat = new ArrayRealVector(singleHyp.Xi);
		    singleHyp.xhatpre = new ArrayRealVector(singleHyp.Xi);
		    singleHyp.Pi = singleHyp.P.copy().getData();
		    singleHyp.weight = 0.5/numsta;
		    singleHyp.spvec = new double[veclen];

		    singleHyp.tm = new AbsoluteDate(singleHyp.epoch, 0.0);
		    singleHyp.ssta = new SpacecraftState[1];
	
		    singleHyp.prop = new ManualPropagation(odcfg, veclen);
		    
		    singleHyp.ProbD = odcfg.Estimation.ProbDetection;
		    singleHyp.smootherResults = singleHyp.new smootherData();
		    singleHyp.hypothesisWeight = 1;
		    
		    ArrayList<SCstate> tempList = new ArrayList<SCstate>();
		    tempList.add(singleHyp);
			states.set(objNum, tempList);

		}
		
		rawMeasurements =  new ArrayList<Measurements.JSONMeasurement>(Arrays.asList(odobs.rawmeas));

		int maxSmootherIter = states.get(0).get(0).odcfg.Estimation.USKFIterations;
		boolean CAREnabled = false;

		for(int smootherIter = 0; smootherIter < maxSmootherIter; smootherIter++)
		{
	    	System.out.printf("Iteration: %d \n", smootherIter);

			unassociatedObsJSON =  new ArrayList<>(rawMeasurements);

			

			
			
			
			if(smootherIter >= 1)
			{
				
				for(int objNum = 0; objNum < states.size(); objNum++)
				{

						states.get(objNum).get(0).xhat = new ArrayRealVector(states.get(objNum).get(0).objResults.Estimation.get(0).EstimatedState);
				    	states.get(objNum).get(0).xhatpre = new ArrayRealVector(states.get(objNum).get(0).objResults.Estimation.get(0).EstimatedState);
				    	states.get(objNum).get(0).P = new DiagonalMatrix(states.get(objNum).get(0).odcfg.Estimation.Covariance);
					    states.get(objNum).get(0).tm = new AbsoluteDate(states.get(objNum).get(0).objResults.Estimation.get(0).Time, DataManager.utcscale);
						states.get(objNum).get(0).smootherResults = states.get(objNum).get(0).new smootherData();
						
						states.get(objNum).get(0).objResults = new JSONResults();
						states.get(objNum).get(0).associatedObsJSON = new ArrayList<Measurements.JSONMeasurement>();
						states.get(objNum).get(0).associatedObs = new ArrayList<ObservedMeasurement>();
						states.get(objNum).get(0).MarginalEvents = new ArrayList<JPDALikelihoods>();
						
				}
			}
			

			
		int mix = -1;
		int additionalMeas = 0;
		
		boolean ODfinished = false;

		//// Loop for all obs starts
		EstimationLoop:
	    while(!ODfinished)
	    {

	    mix = mix + additionalMeas + 1;
    	System.out.printf("Percent Complete: %.2f%% \n", (double) mix/rawMeasurements.size() * 100);

	    //calc additional meas
	    additionalMeas = 0;

		if (mix + additionalMeas < rawMeasurements.size())
		{
		    while(true)
		    {	    	

		    	int i = 1; 
		    	if(mix+additionalMeas+1 < rawMeasurements.size() && rawMeasurements.get(mix+additionalMeas).Time.equals(rawMeasurements.get(mix+additionalMeas+1).Time)  
		    			&& rawMeasurements.get(mix+additionalMeas).Station.equals(rawMeasurements.get(mix+additionalMeas+1).Station))
		    	{
		    		additionalMeas++;
		    		i++;
		    	}
		    	else
		    	{
		    		break;
		    	}
		    }
		}
		else
		{		

			ODfinished = true;
		}
		
		//If CAR is turned on & no sc states currently exist, then create states with earliest measurements
		//Note that since multiple earliest measurements can exist, read in all of the measurements.
		
		if(CAREnabled == true && states.size() == 0)
		{
			for(int measNum = 0; measNum <= additionalMeas; measNum++)
			{
				System.out.println(rawMeasurements.get(measNum));
				System.out.println(promotedTracks.get(0).odcfg.stations.get(rawMeasurements.get(measNum).Station));
				System.out.println(promotedTracks.get(0).odcfg.propframe);
				System.out.println(new AbsoluteDate(DateTimeComponents.parseDateTime(rawMeasurements.get(measNum).Time),DataManager.utcscale));
				System.out.println(Math.sqrt(promotedTracks.get(0).R.getEntry(0,0)));
				System.out.println(promotedTracks.get(0).RaNoise);
				System.out.println(promotedTracks.get(0).DecNoise);
				System.out.println("asd");
				
				ArrayList<SCstate> tempHypotheses = new ArrayList<SCstate>();

				
				ArrayList<Hypothesis> newHypotheses = GenerateHypotheses(rawMeasurements.get(measNum), promotedTracks.get(0).odcfg.stations.get(rawMeasurements.get(measNum).Station), promotedTracks.get(0).odcfg.propframe,
													  new AbsoluteDate(DateTimeComponents.parseDateTime(rawMeasurements.get(measNum).Time),DataManager.utcscale),
													  Math.sqrt(promotedTracks.get(0).R.getEntry(0,0)), promotedTracks.get(0).RaNoise, Math.sqrt(promotedTracks.get(0).R.getEntry(1,1)), promotedTracks.get(0).DecNoise);
				  



				
				for(int hypNum = 0; hypNum < newHypotheses.size(); hypNum++)
				{
					
					SCstate tempState = new SCstate();
					

					tempState.odcfg = promotedTracks.get(0).odcfg;
					tempState.Q = tempState.odcfg.getProcessNoiseMatrix();
					
					tempState.Rsize = 0;
				    for (String s: meanames)
				    	tempState.Rsize += tempState.odcfg.Measurements.get(s).Error.length;

				    tempState.R = new Array2DRowRealMatrix(tempState.Rsize, tempState.Rsize);
				    for (int i = 0, j = 0; i < meanames.length; i++)
				    {
					Settings.JSONMeasurement jm = tempState.odcfg.Measurements.get(meanames[i]);
					for (int k = 0; k < jm.Error.length; k++)
					{
						tempState.R.setEntry(j, j, jm.Error[k]*jm.Error[k]);
					    j++;
					}
				    }

				    tempState.weight = 0.5/numsta;
				    tempState.spvec = new double[veclen];
				    tempState.ssta = new SpacecraftState[1];
				    tempState.prop = new ManualPropagation(tempState.odcfg, veclen);
				    tempState.ProbD = tempState.odcfg.Estimation.ProbDetection;
				    tempState.smootherResults = tempState.new smootherData();
				   					
				    tempState.Xi = tempState.odcfg.getInitialState();
				    
					for(int j=0; j < 6; j++)
					{
						tempState.Xi[j] = newHypotheses.get(hypNum).xhat.getEntry(j);

					}

					tempState.xhat = new ArrayRealVector(tempState.Xi);
					tempState.xhatpre = new ArrayRealVector(tempState.Xi);
					
					
					
					tempState.P = new Array2DRowRealMatrix(tempState.Xi.length, tempState.Xi.length);

					for(int j=0; j < tempState.Xi.length - 6; j++)
					{
							tempState.P.setEntry(6+j, 6+j, odcfg.Estimation.Covariance[6+j]);
				
					}
							
					
					for(int j=0; j < 6; j++)
					{
						for(int k=0; k < 6; k++)
						{
							tempState.P.setEntry(j, k, newHypotheses.get(hypNum).P.getEntry(j,k));

						}					
					}

					tempState.Pi = tempState.P.getData();

					tempState.hypothesisWeight = newHypotheses.get(hypNum).weight;
					
					tempState.epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(rawMeasurements.get(measNum).Time),DataManager.utcscale);				
					tempState.tm = tempState.epoch;
					
					tempState.smootherResults = tempState.new smootherData();
					tempState.objResults = new JSONResults();
					tempState.associatedObsJSON = new ArrayList<Measurements.JSONMeasurement>();
					tempState.associatedObs = new ArrayList<ObservedMeasurement>();
					tempState.dataAssociated = true;

					
					tempHypotheses.add(tempState);

				}
				
				states.add(tempHypotheses);

			}
			continue;
		}
		
		
		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
			{
				states.get(objNum).get(hypNum).MarginalEvents.clear();
			}
		}
		
		ArrayList<ArrayList<JPDALikelihoods>> JointEvents = new ArrayList<ArrayList<JPDALikelihoods>>();
		ArrayList<JPDALikelihoods> SingleJointEvent = new ArrayList<JPDALikelihoods>();

		
		
		///////////////////////////////////////////////////////// FIRST OBJECT LOOP /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////// Propagation & Constructing Marginal Events/////////////////////////////////
		System.out.println("Total Number of Hypotheses : "  + states.get(0).size());

		for(int objNum = 0; objNum < states.size(); objNum++)
		{
		for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
		{
		System.out.println("HypNum:" + hypNum);

		SCstate currentObj = states.get(objNum).get(hypNum);
		odcfg = currentObj.odcfg;
		

		
		currentObj.odout = results.new JSONEstimation();

		Array2DRowRealMatrix sigma = new Array2DRowRealMatrix(numsta, numsig);
		Array2DRowRealMatrix spupd = new Array2DRowRealMatrix(currentObj.Rsize, numsig);
	    
	    
		AbsoluteDate t0 = currentObj.tm;
		
		if(!ODfinished)
		{
			
			currentObj.tm = new AbsoluteDate(DateTimeComponents.parseDateTime(rawMeasurements.get(mix).Time),
					  DataManager.utcscale);

		    
		    double[] pv = currentObj.xhat.toArray();
		    TimeStampedPVCoordinates pvs = new TimeStampedPVCoordinates(currentObj.tm,
										new Vector3D(pv[0], pv[1], pv[2]),
										new Vector3D(pv[3], pv[4], pv[5]));

		    if (odcfg.stations != null && rawMeasurements.get(mix).Station != null)
		    {

			PVCoordinates pvi = odcfg.stations.get(rawMeasurements.get(mix).Station).getBaseFrame().
			    getPVCoordinates(currentObj.tm, odcfg.propframe);
			
		    }

		    states.get(objNum).get(hypNum).odout.Time = rawMeasurements.get(mix).Time;

		}
		else
		{

			currentObj.tm = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.End),
					  DataManager.utcscale);
		    

		}

		
		if(currentObj.dataAssociated == true || mix == 0)
		{
			
		states.get(objNum).get(hypNum).dataAssociated = false;
		
		currentObj.sigpr = new Array2DRowRealMatrix(numsta, numsig);
		
		
		RealMatrix Ptemp = currentObj.P.scalarMultiply(numsta);
		//RealMatrix sqrP = new CholeskyDecomposition(
		//   Ptemp.add(Ptemp.transpose()).scalarMultiply(0.5), 1E-6, 1E-14).getL(); ///////////////////////////////////////////////////////////////////////////////////////////////////////////

		RealMatrix sqrP = new Array2DRowRealMatrix(numsta,numsta);
		try{sqrP = new CholeskyDecomposition(
			    Ptemp.add(Ptemp.transpose()).scalarMultiply(0.5), 1E-6, 1E-14).getL();}
		catch(Exception e)
		{		
			

			System.out.println(e.getStackTrace());
			System.out.println("Ptemp: " + Ptemp);
			System.exit(0);


		}
		
		for (int i = 0; i < numsta; i++)
		{
		    sigma.setColumnVector(i, currentObj.xhat.add(sqrP.getColumnVector(i)));
		    sigma.setColumnVector(numsta + i, currentObj.xhat.subtract(sqrP.getColumnVector(i)));
		}

		if (odcfg.estparams.size() > 0)
		{
		    double[][] sigdata = sigma.getData();

		    for (int j = 6; j < odcfg.estparams.size() + 6; j++)
		    {
			Settings.EstimatedParameter tempep = odcfg.estparams.get(j - 6);
			for (int i = 0; i < numsig; i++)
			    sigdata[j][i] = Math.min(Math.max(sigdata[j][i], tempep.min), tempep.max);
		    }

		    sigma.setSubMatrix(sigdata, 0, 0);
		}
		
		}
		else
		{
			sigma = currentObj.sigpr;

		}


		double propt0 = t0.durationFrom(currentObj.epoch);
		double propt1 = currentObj.tm.durationFrom(currentObj.epoch);

		
		if (propt0 == propt1)
			currentObj.sigpr.setSubMatrix(sigma.getData(), 0, 0);
		else
		{
			if(states.get(objNum).size()>1)
			{
				try{unstack(currentObj.sigpr, currentObj.prop.propagate(propt0, stack(sigma, currentObj.spvec), propt1));}
				catch(OrekitException e)
				{		
					

					currentObj.badHypothesis = true;
					continue;


				}
			}
			else
			    unstack(currentObj.sigpr, currentObj.prop.propagate(propt0, stack(sigma, currentObj.spvec), propt1));
		}

		
		
		
		currentObj.xhatpre = addColumns(currentObj.sigpr).mapMultiplyToSelf(currentObj.weight);

		
		if (ODfinished) 
		{
		    
		    if(objNum == states.size()-1)
		    {
		    	//All objects have been propagated
		    	break EstimationLoop;
		    }
		    else
		    {
		    	//propagate next object to final timestep
		    	continue;
		    }
		    
		    
		}
		
		
		currentObj.Ppre = currentObj.Q.copy();




		for (int i = 0; i < numsig; i++)
		{
		    RealVector y = currentObj.sigpr.getColumnVector(i).subtract(currentObj.xhatpre);
		    currentObj.Ppre = currentObj.Ppre.add(y.outerProduct(y).scalarMultiply(currentObj.weight));

		    double[] pv = currentObj.sigpr.getColumn(i);
		    currentObj.ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])),
								     odcfg.propframe, currentObj.tm, Constants.EGM96_EARTH_MU),
		    		currentObj.prop.getAttitude(currentObj.tm, pv), odcfg.SpaceObject.Mass);


		    
		    if (combmeas)
		    {
				double[] fitv = odobs.measobjs.get(mix).estimate(1, 1, currentObj.ssta).getEstimatedValue();
				spupd.setColumn(i, fitv);
				

		    }
		    else
		    {
				double[] fitv = odobs.measobjs.get(mix*2).estimate(1, 1, currentObj.ssta).getEstimatedValue();
				spupd.setEntry(0, i, fitv[0]);

			if (currentObj.Rsize > 1)
			{
			    fitv = odobs.measobjs.get(mix*2 + 1).estimate(1, 1, currentObj.ssta).getEstimatedValue();
			    spupd.setEntry(1, i, fitv[0]);
			}

		    }
		    
		    
		    
		}

		currentObj.Pyy = currentObj.R.copy();
		RealMatrix Pxy = new Array2DRowRealMatrix(numsta, currentObj.Rsize);
		currentObj.yhatpre = addColumns(spupd).mapMultiplyToSelf(currentObj.weight);
		for (int i = 0; i < numsig; i++)
		{
		    RealVector y = spupd.getColumnVector(i).subtract(currentObj.yhatpre);
		    currentObj.Pyy = currentObj.Pyy.add(y.outerProduct(y).scalarMultiply(currentObj.weight));
		    Pxy = Pxy.add(currentObj.sigpr.getColumnVector(i).subtract(currentObj.xhatpre).outerProduct(y).scalarMultiply(currentObj.weight));
		}


		currentObj.K = Pxy.multiply(MatrixUtils.inverse(currentObj.Pyy));

		//start meas loop


		for(int measNum = 0; measNum <= additionalMeas; measNum++)
		{
			
			RealVector raw = null;
			//adjust code to recompute raw	

			    if (combmeas)
			    {
					    raw = new ArrayRealVector(odobs.measobjs.get(mix+measNum).getObservedValue());
			    }
			    else
			    {

					if (currentObj.Rsize > 1)
					{
						raw = new ArrayRealVector(new double[]{odobs.measobjs.get((mix+measNum)*2).getObservedValue()[0],
										       odobs.measobjs.get((mix+measNum)*2 + 1).getObservedValue()[0]});
					}
					else
					{
					    raw = new ArrayRealVector(new double[]{odobs.measobjs.get((mix+measNum)*2).getObservedValue()[0]});
					}
			    }

			   
			// add "measurement 0" cases to jpda object
			if(currentObj.MarginalEvents.isEmpty())
			{
				JPDALikelihoods JPDAtemp = new JPDALikelihoods();
				
				JPDAtemp.Psi = currentObj.hypothesisWeight * (1 - currentObj.ProbD);
				JPDAtemp.object = objNum;
				JPDAtemp.measurement = 0;

				states.get(objNum).get(hypNum).MarginalEvents.add(JPDAtemp);


				
				
			}
		


		// run through jpda algorithms and calculate psi_r_i
			
	    //Innovations QuadCheck
			
		RealVector Innov = quadCheck(raw, currentObj.yhatpre);
	

							
		RealMatrix MahalaTemp = MatrixUtils.createColumnRealMatrix(Innov.toArray());

		RealMatrix Mahalanobis = MahalaTemp.transpose().multiply(MatrixUtils.inverse(currentObj.Pyy).multiply(MahalaTemp));


		
		
		if( odcfg.Estimation.GatingThreshold > Math.sqrt(Mahalanobis.getEntry(0,0)))
		{
						
			JPDALikelihoods JPDAtemp = new JPDALikelihoods();

			//JPDAtemp2.Psi =  Math.exp(-JPDATemp.getEntry(0,0)/2) / (Math.sqrt(new LUDecomposition(currentObj.Pyy).getDeterminant() * Math.pow(2 * Math.PI, currentObj.Rsize)));
			JPDAtemp.Psi = currentObj.hypothesisWeight * (1 - new ChiSquaredDistribution(currentObj.Rsize).cumulativeProbability(Mahalanobis.getEntry(0,0)));
						
			JPDAtemp.object = objNum;
			JPDAtemp.measurement = measNum+1;

			
			currentObj.MarginalEvents.add(JPDAtemp);
			

			
		}

					
		// end loops (measurement then hypothesis then object)
		}
		}

	    }
		
		
		
		//remove BadHypotheses and normalize.
		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			int hypNum = 0;
			
			double sumWeights = 0;
			
			while(hypNum < states.get(objNum).size())
			{	
				if(states.get(objNum).get(hypNum).badHypothesis == true)
				{
					states.get(objNum).remove(hypNum);
				}
				else
				{
					sumWeights += states.get(objNum).get(hypNum).hypothesisWeight;
					hypNum++;
				}
			}
			
			for(hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
			{
				states.get(objNum).get(hypNum).hypothesisWeight /= sumWeights;
				
				for(int eventCounter = 0; eventCounter < states.get(objNum).get(hypNum).MarginalEvents.size(); eventCounter++)
				{
					states.get(objNum).get(hypNum).MarginalEvents.get(eventCounter).Psi /= sumWeights;
				}
			}
		}
		
		
		/////////////////////////////////////////////////////////////////////////////////////////
		if(CAREnabled)
		{
		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			
			try {
				FileWriter writer = new FileWriter("../../../Matlab/outputStates.txt",true); 

				for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
				{
					SCstate tempState = states.get(objNum).get(hypNum);
					for(int i = 0; i <  tempState.xhatpre.getDimension(); i++)
					{
						writer.write(tempState.xhatpre.getEntry(i)+",");
					}
					
					for(int i = 0; i <  tempState.Ppre.getRowDimension(); i++)
					{
						for(int j = 0; j <  tempState.Ppre.getColumnDimension(); j++)
						{
							writer.write(tempState.Ppre.getEntry(i,j)+",");
						}
					}
				
				
			  writer.write(tempState.hypothesisWeight+"," +mix + System.lineSeparator());
			  }
			
			writer.close();
		
			} catch (IOException e) {
		
			e.printStackTrace();
		
			} 
		}
		}
		/////////////////////////////////////////////////////////////////////////////////////////
		

		
		//Determine joint events
		double[][] SumJPDALikelihoods = new double[states.size()][additionalMeas+2]; // +1 to account for measurement 0 case
			
		
		for(int objNum = 0; objNum < states.size(); objNum++)
		{

			for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
			{
				SCstate currentHyp = states.get(objNum).get(hypNum);
				
				
				for(int eventCounter = 0; eventCounter < currentHyp.MarginalEvents.size(); eventCounter++)
				{
					
					SumJPDALikelihoods[objNum][currentHyp.MarginalEvents.get(eventCounter).measurement] += currentHyp.MarginalEvents.get(eventCounter).Psi;
										
				}
				
			}		
			
		}

///////////////////////////////////////////////////
		if(CAREnabled)
		{
			for(int objNum = 0; objNum < states.size(); objNum++)
			{
			
				try {
					
					FileWriter writer = new FileWriter("../../../Matlab/MarginalEventPsi.txt",true); 						
						
						for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
						{
							SCstate currentHyp = states.get(objNum).get(hypNum);

							for(int eventCounter = 0; eventCounter < currentHyp.MarginalEvents.size(); eventCounter++)
							{
								writer.write(mix+","+objNum+","+hypNum+","+currentHyp.MarginalEvents.get(eventCounter).measurement+","+currentHyp.MarginalEvents.get(eventCounter).Psi);
								writer.write(System.lineSeparator());
							}
						}
						
					
						writer.close();
				
				} catch (IOException e) {
				
					e.printStackTrace();
				
				} 
			}
		}	
		
///////////////////////////////////////////////////
	if(CAREnabled)
	{
		for(int objNum = 0; objNum < SumJPDALikelihoods.length; objNum++)
		{
		
			try {
				
				FileWriter writer = new FileWriter("../../../Matlab/SumJPDALiklihood.txt",true); 
	
					writer.write(mix+","+objNum+","+SumJPDALikelihoods[objNum].length+",");
					
					
					for(int measNum = 0; measNum < SumJPDALikelihoods[objNum].length; measNum++)
					{
						writer.write(SumJPDALikelihoods[objNum][measNum]+",");
					}
					
					writer.write(System.lineSeparator());
				
				
					writer.close();
			
			} catch (IOException e) {
			
				e.printStackTrace();
			
			} 
		}
	}
///////////////////////////////////////////////////
		

		JointEvents = JPDAJointEvents(JointEvents, SumJPDALikelihoods, SingleJointEvent, 0);
		
		//compute Probabilities
		double[] JPDAProbability = new double[JointEvents.size()];
		Arrays.fill(JPDAProbability, 1);
		
		for(int i = 0; i<JointEvents.size(); i++)
		{
			for(int j = 0; j<JointEvents.get(i).size(); j++)
			{
				//skip cases where  only one marginal event to choose from. This is important for when Pd = 1, and an object does not 
				//have an associated observation. This is detected by checking if an object has a probability of association of 0 for all measurements (including 0).
				//If this test is not checked then all joint probabilities will be 0
				
				double rowSum = 0;
				
				for(int measNum = 0; measNum < SumJPDALikelihoods[j].length; measNum++)
				{
					rowSum += SumJPDALikelihoods[j][measNum];
				}
				
				
				if(rowSum > 0)
				{
					JPDAProbability[i] = JPDAProbability[i] * JointEvents.get(i).get(j).Psi;
				}
				
						
			}	
		}


		
		double JPDAsum= Arrays.stream(JPDAProbability).sum();
				
				
		for(int i = 0; i<JointEvents.size(); i++)
		{
			JPDAProbability[i] = JPDAProbability[i] / JPDAsum;
		}	

		
		//identify max probability to see if object has been associated
		int maxProbIndex = 0;
		
		for(int i=1;i < JPDAProbability.length;i++)
		{
		    if(JPDAProbability[i] > JPDAProbability[maxProbIndex])
		    	maxProbIndex = i;
		}


		for(int i = 0; i < JointEvents.get(maxProbIndex).size(); i++)
		{	
			if(JointEvents.get(maxProbIndex).get(i).measurement != 0)
			{
				//remove from unassociated
				unassociatedObsJSON.remove(rawMeasurements.get(mix + JointEvents.get(maxProbIndex).get(i).measurement - 1));


				if(states.get(JointEvents.get(maxProbIndex).get(i).object).size() == 1)
				{
					//output
					states.get(JointEvents.get(maxProbIndex).get(i).object).get(0).associatedObsJSON.add( 
							rawMeasurements.get(mix + JointEvents.get(maxProbIndex).get(i).measurement - 1));
					
					//used for Innov smoother later
					if(combmeas)
					{
						states.get(JointEvents.get(maxProbIndex).get(i).object).get(0).associatedObs.add( 
								odobs.measobjs.get(mix + JointEvents.get(maxProbIndex).get(i).measurement - 1));			
					}
					else
					{
						states.get(JointEvents.get(maxProbIndex).get(i).object).get(0).associatedObs.add( 
								odobs.measobjs.get(2 * (mix + JointEvents.get(maxProbIndex).get(i).measurement - 1)));
						
						states.get(JointEvents.get(maxProbIndex).get(i).object).get(0).associatedObs.add( 
								odobs.measobjs.get(2 * (mix + JointEvents.get(maxProbIndex).get(i).measurement - 1) + 1));					
						
		
					}
				}
				
				for(int j = 0; j < states.get(JointEvents.get(maxProbIndex).get(i).object).size(); j++)
				states.get(JointEvents.get(maxProbIndex).get(i).object).get(j).dataAssociated = true; //indicate all hypotheses for a given object have been associated.
			}
		
		}


		double[][] beta_sat_meas = new double[states.size()][additionalMeas+2]; // +1 to account for measurement 0 case
		
		///////////////////////////////////////////////////////// SECOND OBJECT LOOP /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////// Combine hypotheses into events based upon object & measurement ///////////////////

		
		
		for(int objNum = 0; objNum < states.size(); objNum++)
		{	
			
			SCstate currentObj = states.get(objNum).get(0);
			odcfg = currentObj.odcfg;

			for(int measNum = 0; measNum <=additionalMeas; measNum++)
			{
			
				//compute Beta based on the object & measurement pair

				for(int i = 0; i < JointEvents.size(); i++)
				{

					for(int j = 0; j < JointEvents.get(i).size(); j++)
					{

						if(JointEvents.get(i).get(j).object == objNum && JointEvents.get(i).get(j).measurement == 0)
						{
							
						    beta_sat_meas[objNum][0] += JPDAProbability[i];

						}
						
						if(JointEvents.get(i).get(j).object == objNum && JointEvents.get(i).get(j).measurement - 1 == measNum)
						{
							
						    beta_sat_meas[objNum][measNum+1] += JPDAProbability[i];

						}
					}	
				}		
			}	
		}

		
///////////////////////////////////////////////////
		if(CAREnabled)
		{
		for(int objNum = 0; objNum < beta_sat_meas.length; objNum++)
		{

				try {
					
					FileWriter writer = new FileWriter("../../../Matlab/betaSatMeas.txt",true); 

					
						writer.write(mix+","+objNum+","+beta_sat_meas[objNum].length+",");
	
		
						for(int measNum = 0; measNum < beta_sat_meas[objNum].length; measNum++)
						{
							writer.write(beta_sat_meas[objNum][measNum]+",");
						}
						
						  writer.write(System.lineSeparator());

				
				writer.close();
			
				} catch (IOException e) {
			
				e.printStackTrace();
			
				} 
		}
		}
///////////////////////////////////////////////////

		
		///////////////////////////////////////////////////////// THIRD OBJECT LOOP /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////// Compute probability of joint events, Perform update step ///////////////////		

		for(int objNum = 0; objNum < states.size(); objNum++)
		{

		for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
		{	

			SCstate currentObj = states.get(objNum).get(hypNum);
			odcfg = currentObj.odcfg;
			
			double Beta = 0;
			RealMatrix Betav = new Array2DRowRealMatrix(currentObj.Rsize, 1);
			RealMatrix Betavvt = new Array2DRowRealMatrix(currentObj.Rsize, currentObj.Rsize);

			//update hypothesis weight
			if(states.get(objNum).size() > 1)
			{
				double hypWeight = 0;			
				
				for(int eventCounter = 0; eventCounter < currentObj.MarginalEvents.size(); eventCounter++)
				{
					
					int measNum = currentObj.MarginalEvents.get(eventCounter).measurement;
					
					hypWeight += beta_sat_meas[objNum][measNum] / SumJPDALikelihoods[objNum][measNum] * currentObj.MarginalEvents.get(eventCounter).Psi;
				}
				
				states.get(objNum).get(hypNum).hypothesisWeight = hypWeight;
				
///////////////////////////////////////////////////
				if(CAREnabled)
				{

					try {
						
						FileWriter writer = new FileWriter("../../../Matlab/HypWeight.txt",true); 
			
							writer.write(mix+","+objNum+","+hypNum+","+states.get(objNum).size()+","+states.get(objNum).get(hypNum).hypothesisWeight);
							
							writer.write(System.lineSeparator());
						
							writer.close();
					
					} catch (IOException e) {
					
						e.printStackTrace();
					
					} 
				}
///////////////////////////////////////////////////
			}

			
			for(int eventCounter = 0; eventCounter < currentObj.MarginalEvents.size(); eventCounter++)
			{
			
				//compute Beta
				//check each probability to see if it has the corresponding obj/meas combo

				if(currentObj.MarginalEvents.get(eventCounter).measurement > 0)
				{
					int measNum = currentObj.MarginalEvents.get(eventCounter).measurement - 1;
					RealVector raw = null;

				    if (combmeas)
				    {
						    raw = new ArrayRealVector(odobs.measobjs.get(mix+measNum).getObservedValue());
				    }
				    else
				    {

						if (currentObj.Rsize > 1)
						{
							raw = new ArrayRealVector(new double[]{odobs.measobjs.get((mix+measNum)*2).getObservedValue()[0],
											       odobs.measobjs.get((mix+measNum)*2 + 1).getObservedValue()[0]});
						}
						else
						{
						    raw = new ArrayRealVector(new double[]{odobs.measobjs.get((mix+measNum)*2).getObservedValue()[0]});
						}
				    }
						    
						    
				    //Innovations QuadCheck
				    RealMatrix Innovation = MatrixUtils.createColumnRealMatrix(quadCheck(raw, currentObj.yhatpre).toArray());

				    double BetaTemp = 0;
				    
				    //This if statement is necessary if sumJPDALikelihoods is 0. This can only happen when JointEvents.get(i).get(j).Psi ==0 for all j. Betatemp now defaults to
				    //0, and if JointEvents.get(i).get(j).Psi != 0, then it updates. If JointEvents.get(i).get(j).Psi == 0, then Betatemp should be 0 regardless
				    //This check avoids a 0/0 case.
				    if(currentObj.MarginalEvents.get(eventCounter).Psi != 0)
				    {

					    BetaTemp = beta_sat_meas[objNum][measNum+1] / SumJPDALikelihoods[objNum][measNum+1] * currentObj.MarginalEvents.get(eventCounter).Psi / currentObj.hypothesisWeight;
					    

				    }

				    


					Beta = Beta + BetaTemp;
					Betav= Betav.add(Innovation.scalarMultiply(BetaTemp));
					Betavvt = Betavvt.add(Innovation.multiply(Innovation.transpose()).scalarMultiply(BetaTemp));


				}	
				
			}	
			

			
			
			// update step

			int numConsideredParams = odcfg.considerparams.size();
			int numEstimatedParams = numsta - numConsideredParams;

			RealMatrix KnoConsider = new Array2DRowRealMatrix(numsta, currentObj.Rsize);

			for(int i = 0; i<numEstimatedParams; i++)
			{	 		

				KnoConsider.setRowVector(i,currentObj.K.getRowVector(i));
			}


			currentObj.xhat = new ArrayRealVector(currentObj.xhatpre.add(KnoConsider.multiply(Betav).getColumnVector(0)));

			
			for (int i = 0; i < (numsta-numEstimatedParams); i++)
		    {
		        currentObj.xhat.setEntry(i+numEstimatedParams, currentObj.Xi[i+numEstimatedParams]);
		    }

			RealMatrix Ptilda = currentObj.K.multiply((Betavvt.subtract(Betav.multiply(Betav.transpose()))).multiply(currentObj.K.transpose()));
		    RealMatrix Pcorrective = currentObj.K.multiply(currentObj.Pyy.multiply(currentObj.K.transpose())).scalarMultiply(Beta).subtract(Ptilda);    //Ptilda is subtracted because later PNoconsider is subtracted. Ptilda         
		    RealMatrix PcorrectiveNoConsider = new Array2DRowRealMatrix(numsta,numsta);                                // needs to be included in the no consider part to remove update to considered params
	


			
		    for (int i = 0; i < numEstimatedParams; i++)          
		    {
		        PcorrectiveNoConsider.setRowVector(i, Pcorrective.getRowVector(i));
		        PcorrectiveNoConsider.setColumnVector(i, Pcorrective.getColumnVector(i));
		    }
	
			
			currentObj.P = currentObj.Ppre.subtract(PcorrectiveNoConsider);
			
			

			double[] pv = currentObj.xhat.toArray();

			
			currentObj.ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
											   new Vector3D(pv[3], pv[4], pv[5])),
									 odcfg.propframe, currentObj.tm, Constants.EGM96_EARTH_MU),
						      currentObj.prop.getAttitude(currentObj.tm, pv), odcfg.SpaceObject.Mass);
			if (currentObj.ssta[0].getA() <= Constants.WGS84_EARTH_EQUATORIAL_RADIUS)
			    throw(new RuntimeException(String.format("Invalid semi-major axis %f", currentObj.ssta[0].getA())));
	
			currentObj.odout.EstimatedState = pv;
			currentObj.odout.EstimatedCovariance = currentObj.P.getData();
			currentObj.odout.InnovationCovariance = currentObj.Pyy.getData();
			currentObj.odout.numConsideredParams = numConsideredParams;


						

			
			if (combmeas)
			{
			    for (int i = 0; i < meanames.length; i++)
			    {	
				double[] fitv = odobs.measobjs.get(mix).estimate(1, 1, currentObj.ssta).getEstimatedValue();
				
				if (meanames.length == 1)
				{
					currentObj.odout.PreFit.put(meanames[i], currentObj.yhatpre.toArray());
					currentObj.odout.PostFit.put(meanames[i], fitv);
				}
				else
				{
					currentObj.odout.PreFit.put(meanames[i], new double[] {currentObj.yhatpre.getEntry(i)});
					currentObj.odout.PostFit.put(meanames[i], new double[] {fitv[i]});
				    
				    
				}
			    }
			    
			    
			}
			else
			{

			    double[] fitv = odobs.measobjs.get(mix*2).estimate(1, 1, currentObj.ssta).getEstimatedValue();
			    currentObj.odout.PreFit.put(meanames[0], new double[] {currentObj.yhatpre.getEntry(0)});
			    currentObj.odout.PostFit.put(meanames[0], fitv);
	
			    if (currentObj.Rsize > 1)
			    {
				fitv = odobs.measobjs.get(mix*2 + 1).estimate(1, 1, currentObj.ssta).getEstimatedValue();
				currentObj.odout.PreFit.put(meanames[1], new double[] {currentObj.yhatpre.getEntry(1)});
				currentObj.odout.PostFit.put(meanames[1], fitv);
			    }
			}
			
			if(currentObj.dataAssociated == true)
			currentObj.objResults.Estimation.add(currentObj.odout);


			
			
			if(odcfg.Estimation.Smoother.equals("On") && currentObj.dataAssociated == true)
			{	
				
				
			// compute post sigma points
		    Array2DRowRealMatrix sigpost = new Array2DRowRealMatrix(numsta, numsig);
	
			RealMatrix Ptemp2 = currentObj.P.scalarMultiply(numsta);
		    
			
			RealMatrix sqrP2 = new CholeskyDecomposition(
			    Ptemp2.add(Ptemp2.transpose()).scalarMultiply(0.5), 1E-6, 1E-14).getL();

			//generate points
			for (int i = 0; i < numsta; i++)
			{
				sigpost.setColumnVector(i, currentObj.xhat.add(sqrP2.getColumnVector(i)));
				sigpost.setColumnVector(numsta + i, currentObj.xhat.subtract(sqrP2.getColumnVector(i)));
			}
			
			//reduce estimated params to max/min
			if (odcfg.estparams.size() > 0)
			{
			    double[][] sigdata = sigpost.getData();
	
			    for (int j = 6; j < odcfg.estparams.size() + 6; j++)
			    {
				Settings.EstimatedParameter tempep = odcfg.estparams.get(j - 6);
				for (int i = 0; i < numsig; i++)
				    sigdata[j][i] = Math.min(Math.max(sigdata[j][i], tempep.min), tempep.max);
			    }
	
			    sigpost.setSubMatrix(sigdata, 0, 0);
			}
			
			// store smoother data

				SCstate.smootherData.smootherStep smout = currentObj.smootherResults.new smootherStep();
		
				smout.xpre = MatrixUtils.createColumnRealMatrix(currentObj.xhatpre.toArray());
				smout.xpost = MatrixUtils.createColumnRealMatrix(currentObj.xhat.toArray());
				
				smout.Ppre = currentObj.Ppre;
				smout.Ppost = currentObj.P;
		
				smout.sigPre = MatrixUtils.createRealMatrix(currentObj.sigpr.getData());
				smout.sigPost = MatrixUtils.createRealMatrix(sigpost.getData());
		
				smout.tmSmoother = currentObj.tm;
	    		
	    		
				if(combmeas)
				{
					smout.measObjsSmoother = odobs.measobjs.get(mix);
					
				}
				else
				{
					smout.measObjsSmoother = odobs.measobjs.get(2 * mix);

					smout.measObjsSmootherNoComb = odobs.measobjs.get(2* mix + 1);
				}
	    		
				
	    		currentObj.smootherResults.smoother.add(smout);
				
			
			}
				
		    

			
		//end hyp loop then obj
		}
		}
		
		//Remove hypotheses with low weight.
		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			double maxWeight = 0;
			for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
			{
				if(states.get(objNum).get(hypNum).hypothesisWeight > maxWeight)
					maxWeight = states.get(objNum).get(hypNum).hypothesisWeight;
			
			}
			
			
		int hypNum = 0;
				
		while(hypNum < states.get(objNum).size())
		{	
			if(states.get(objNum).get(hypNum).hypothesisWeight < maxWeight * 1E-3)
			{
				states.get(objNum).remove(hypNum);
			}
			else
			{
				hypNum++;
			}
		}
		
		}

		//Merge Hypotheses with low mahalanobis distance
		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			
		int hypNum1 = 0;
				
		while(hypNum1 < states.get(objNum).size())
		{	
			int hypNum2 = hypNum1+1;

			while(hypNum2 < states.get(objNum).size())
			{	
								
				SCstate obj1 = states.get(objNum).get(hypNum1);
				SCstate obj2 = states.get(objNum).get(hypNum2);
				//Combine distributions
				RealMatrix MahalaTemp = MatrixUtils.createColumnRealMatrix(obj1.xhat.subtract(obj2.xhat).toArray());

				RealMatrix Mahalanobis1 = MahalaTemp.transpose().multiply(MatrixUtils.inverse(obj1.P).multiply(MahalaTemp));
				RealMatrix Mahalanobis2 = MahalaTemp.transpose().multiply(MatrixUtils.inverse(obj2.P).multiply(MahalaTemp));
								
				//If mahalanobis distance low, combine distributions.
				if(Math.min(Mahalanobis1.getEntry(0,0),Mahalanobis2.getEntry(0,0)) < 1)
				{
					RealVector xhatTemp =  obj1.xhat.mapMultiply(obj1.hypothesisWeight).add(obj2.xhat.mapMultiply(obj2.hypothesisWeight));
					
					RealMatrix innovTemp = MatrixUtils.createColumnRealMatrix(obj1.xhat.subtract(xhatTemp).toArray());
					
					RealMatrix PTemp = (obj1.P.add(innovTemp.multiply(innovTemp.transpose()))).scalarMultiply(obj1.hypothesisWeight);
					
					innovTemp = MatrixUtils.createColumnRealMatrix(obj2.xhat.subtract(xhatTemp).toArray());

					PTemp = PTemp.add((obj2.P.add(innovTemp.multiply(innovTemp.transpose()))).scalarMultiply(obj2.hypothesisWeight));

					states.get(objNum).get(hypNum1).xhat = new ArrayRealVector(xhatTemp);
					states.get(objNum).get(hypNum1).P = PTemp;
					states.get(objNum).get(hypNum1).hypothesisWeight += states.get(objNum).get(hypNum2).hypothesisWeight;
					
					states.get(objNum).remove(hypNum2);

				}
				else
				{
					hypNum2++;
				}
			}
			
			hypNum1++;
		}
		
		}

		//normalize weights
		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			double sumWeights = 0;
			
			for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
			{	
				sumWeights += states.get(objNum).get(hypNum).hypothesisWeight;
			}
			
			for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
			{
				states.get(objNum).get(hypNum).hypothesisWeight /= sumWeights;
				
			}
		}
		
	    
	    }//while
		
		///////////////////////////////////////////////////////// FOURTH OBJECT LOOP /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////// Smooth Data ///////////////////////////////////////////////////////////////

		for(int objNum = 0; objNum < states.size(); objNum++)
		{

		for(int hypNum = 0; hypNum < states.get(objNum).size(); hypNum++)
		{
			
			SCstate currentObj = states.get(objNum).get(hypNum);
			odcfg = currentObj.odcfg;
			
		    double[] pv = currentObj.xhatpre.toArray();


		    currentObj.tm = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.End),
					  DataManager.utcscale);

		    CartesianOrbit cart  = new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										new Vector3D(pv[3], pv[4], pv[5])),
							      odcfg.propframe, currentObj.tm, Constants.EGM96_EARTH_MU);
		    if (cart.getA() <= Constants.WGS84_EARTH_EQUATORIAL_RADIUS)
			throw(new RuntimeException(String.format("Invalid semi-major axis %f", cart.getA())));

		    
		    
		    currentObj.objResults.Propagation.Time = odcfg.Propagation.End;
		    currentObj.objResults.Propagation.State = pv;
		    

		    
		    // run smoother
				
			currentObj.McReynoldsConsistencyPass = true;
			
			int smSize = currentObj.smootherResults.smoother.size()-1;

			currentObj.smootherResults.smoother.get(smSize).xstar 
				= currentObj.smootherResults.smoother.get(smSize).xpost;
			currentObj.smootherResults.smoother.get(smSize).Pstar 
				= currentObj.smootherResults.smoother.get(smSize).Ppost;			
			

		    for(int i = 0; i < smSize; i++)
		    {
		    	SCstate.smootherData.smootherStep smDatak1 = currentObj.smootherResults.smoother.get(smSize - i);
		    	SCstate.smootherData.smootherStep smDatak = currentObj.smootherResults.smoother.get(smSize - i - 1);
	
				RealMatrix Csmoother = new Array2DRowRealMatrix(numsta,numsta);
				RealMatrix Asmoother = new Array2DRowRealMatrix(numsta,numsta);
		    	
				for(int j = 0; j < numsig; j++) 
				{
					Csmoother = Csmoother.add(((smDatak.sigPost.getColumnMatrix(j).subtract(smDatak.xpost)).multiply( 
					    (smDatak1.sigPre.getColumnMatrix(j).subtract(smDatak1.xpre)).transpose())).scalarMultiply(currentObj.weight));
	
	
				}
		    		
				
				Asmoother = Csmoother.multiply(MatrixUtils.inverse(smDatak1.Ppre));
	
				smDatak.xstar = smDatak.xpost.add(Asmoother.multiply(smDatak1.xstar.subtract(smDatak1.xpre)));
				smDatak.Pstar = smDatak.Ppost.add(Asmoother.multiply((smDatak1.Pstar.subtract(smDatak1.Ppre)).multiply(Asmoother.transpose())));
				


				currentObj.smootherResults.smoother.set(smSize - i - 1,smDatak);
				

		    }
		    
		}

		// compute McReynolds Consistency or merge all hypotheses into one hypothesis
		
		if(states.get(objNum).size() == 1)
		{
			
			double McReynoldsVal = -99999;

			SCstate currentObj = states.get(objNum).get(0);
			int smSize = currentObj.smootherResults.smoother.size()-1;
			
		    for(int i = 0; i < currentObj.smootherResults.smoother.size()-1;i++)
		    {

		    	SCstate.smootherData.smootherStep smDatak1 = currentObj.smootherResults.smoother.get(smSize - i);
		    	SCstate.smootherData.smootherStep smDatak = currentObj.smootherResults.smoother.get(smSize - i - 1);
	
				
				RealMatrix delx = smDatak.xpost.subtract(smDatak.xstar);
						
				RealMatrix delP = smDatak.Ppost.subtract(smDatak.Pstar);
				
				for(int j = 0; j < 5; j++)
				{			

					McReynoldsVal = Math.max(McReynoldsVal, Math.abs(delx.getEntry(j,0)) / Math.sqrt(Math.abs(delP.getEntry(j,j))));
					
					
					if(Math.abs(delx.getEntry(j,0)) / Math.sqrt(Math.abs(delP.getEntry(j,j))) >= 3)
					{


						currentObj.McReynoldsConsistencyPass = false;

					}
				}

			}
			
			
		    System.out.println("McReynolds Pass: " + states.get(objNum).get(0).McReynoldsConsistencyPass +  "   McReynoldsVal: " + McReynoldsVal);
    		System.out.println("Object " + objNum + " Observations Associated: " + states.get(objNum).get(0).associatedObsJSON.size());
    		
    		
    		//store in results
		    
		    SpacecraftState[] smSsta = new SpacecraftState[1];		
		    
		    for(int j = 0; j < currentObj.smootherResults.smoother.size();j++)
		    {
		    	
		    	double[] pv = currentObj.smootherResults.smoother.get(j).xstar.getColumn(0);
			    
			    currentObj.tm = currentObj.smootherResults.smoother.get(j).tmSmoother;
			    
			    smSsta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
											       new Vector3D(pv[3], pv[4], pv[5])), odcfg.propframe, currentObj.tm, Constants.EGM96_EARTH_MU));
	
			    
		    	//compute smoothed residuals
			    
			    if (combmeas)
				{
				    for (int i = 0; i < meanames.length; i++)
				    {
					double[] fitv = currentObj.smootherResults.smoother.get(j).measObjsSmoother.estimate(1, 1, smSsta).getEstimatedValue();
					if (meanames.length == 1)
					{

						currentObj.objResults.Estimation.get(j).PostFit.put(meanames[i], fitv);
					}
					else
					{

						currentObj.objResults.Estimation.get(j).PostFit.put(meanames[i], new double[] {fitv[i]});
						
	
					}
				    }
				}
				else
				{
				    double[] fitv = currentObj.smootherResults.smoother.get(j).measObjsSmoother.estimate(1, 1, smSsta).getEstimatedValue();
				    currentObj.objResults.Estimation.get(j).PostFit.put(meanames[0], fitv);
	
				    if (currentObj.Rsize > 1)
				    {

					fitv = currentObj.smootherResults.smoother.get(j).measObjsSmootherNoComb.estimate(1, 1, smSsta).getEstimatedValue();
					currentObj.objResults.Estimation.get(j).PostFit.put(meanames[1], fitv);
				    }
				}

				
				
		    	//store
			    currentObj.objResults.Estimation.get(j).EstimatedState = currentObj.smootherResults.smoother.get(j).xstar.getColumn(0);
			    currentObj.objResults.Estimation.get(j).EstimatedCovariance = currentObj.smootherResults.smoother.get(j).Pstar.getData();
		    		    	
		    }
    		
		}
		else
		{
    		RealMatrix xstarCombined = states.get(objNum).get(0).smootherResults.smoother.get(0).xstar.copy();
    		RealMatrix PstarCombined = new Array2DRowRealMatrix();
			
			for(int hypNum = 1; hypNum < states.size(); hypNum++)
			{
				SCstate currentObj = states.get(objNum).get(hypNum);
				
				xstarCombined.add(currentObj.smootherResults.smoother.get(0).xstar.scalarMultiply(currentObj.hypothesisWeight));

			}
			
			for(int hypNum = 1; hypNum < states.size(); hypNum++)
			{
				SCstate currentObj = states.get(objNum).get(hypNum);
				
				RealMatrix PstarTemp = currentObj.smootherResults.smoother.get(0).xstar.subtract(xstarCombined);
				
				PstarCombined.add(currentObj.smootherResults.smoother.get(0).Pstar.add(PstarTemp.multiply(PstarTemp.transpose())).scalarMultiply(currentObj.hypothesisWeight));

			}
			
			states.get(objNum).get(0).objResults.Estimation.get(0).EstimatedState = xstarCombined.copy().getColumn(0);
			states.get(objNum).get(0).objResults.Estimation.get(0).EstimatedCovariance = PstarCombined.copy().getData();

			
			states.get(objNum).get(0).Xi = xstarCombined.copy().getColumn(0);
			states.get(objNum).get(0).Pi = PstarCombined.copy().getData();
			
		    System.out.println("Object " + objNum + ": " + states.get(objNum).size() + " hypotheses converged into a single hypothesis");
				
		}

		}

			//break out of smoother iteration if all data has been established
			boolean promotedTrackReset= false;
			int objNum = 0;
			while(objNum < states.size())
			{
				
				if(smootherIter >= 1 && states.get(objNum).size() == 1 && states.get(objNum).get(0).McReynoldsConsistencyPass == true)
				{
					promotedTrackReset = true;
					
					//promote track
					promotedTracks.add(states.get(objNum).get(0));
					
					//remove associated obs from obs dataset
					for(int measObj = 0; measObj < states.get(objNum).get(0).associatedObs.size(); measObj++)
					{
						rawMeasurements.remove(odobs.measobjs.indexOf(states.get(objNum).get(0).associatedObs.get(measObj)));

						odobs.measobjs.remove(states.get(objNum).get(0).associatedObs.get(measObj));
					}
					
					
					states.remove(states.get(objNum));
					
				}
				else
				{
					objNum ++;
				}

			}
			
			if(promotedTrackReset == true)
			{
				if(rawMeasurements.size() == 0)
				{

					break;

				}
				else
				{
					smootherIter = -1; // reset the smoother iteration counter, since some data has been removed.

					for(objNum = 0; objNum < states.size(); objNum++)
					{
						
						states.get(objNum).get(0).xhat = new ArrayRealVector(states.get(objNum).get(0).Xi);
				    	states.get(objNum).get(0).xhatpre = new ArrayRealVector(states.get(objNum).get(0).Xi);
				    	states.get(objNum).get(0).P = new Array2DRowRealMatrix(states.get(objNum).get(0).Pi);
					    states.get(objNum).get(0).tm = new AbsoluteDate(states.get(objNum).get(0).epoch, 0.0);
					    
						states.get(objNum).get(0).smootherResults = states.get(objNum).get(0).new smootherData();
						states.get(objNum).get(0).objResults = new JSONResults();
						states.get(objNum).get(0).associatedObsJSON = new ArrayList<Measurements.JSONMeasurement>();
						states.get(objNum).get(0).associatedObs = new ArrayList<ObservedMeasurement>();
						
						for(int hypNum = 1; hypNum < states.get(objNum).size(); hypNum++)
						{
							states.get(objNum).remove(hypNum);
						}
					}
					
				}
			}

			
			if(CAREnabled == false && (promotedTrackReset == true || states.size() == 0 || smootherIter == maxSmootherIter))
			{
				if(states.size() == 0 || smootherIter == maxSmootherIter)
				CAREnabled = true;
				
				smootherIter = -1;
				continue;
			}

		}

		for(int objNum = 0; objNum < promotedTracks.size(); objNum++)
		{
			SCstate currentObj = promotedTracks.get(objNum);
		    //compute smoothed innov residuals && Quad Check
			
			
			for(int measNum = 0; measNum < currentObj.associatedObs.size(); measNum ++)
			{
				
				
				

				RealMatrix Ptemp = new Array2DRowRealMatrix(currentObj.objResults.Estimation.get(measNum).EstimatedCovariance).scalarMultiply(numsta);
				RealVector xhatTemp = new ArrayRealVector(currentObj.objResults.Estimation.get(measNum).EstimatedState);
				
				RealMatrix sigmaTemp = GenerateSigmaPoints(xhatTemp, Ptemp, numsta, numsig);


				Array2DRowRealMatrix spupd = new Array2DRowRealMatrix(promotedTracks.get(0).Rsize, numsig);

				
				AbsoluteDate tmTemp = currentObj.associatedObs.get(measNum).getDate();

				for (int i = 0; i < numsig; i++)
				{

				    double[] pv = sigmaTemp.getColumn(i);
				    
				    SpacecraftState[] sstaTemp = new SpacecraftState[1];
				    sstaTemp[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
												       new Vector3D(pv[3], pv[4], pv[5])),
										     odcfg.propframe, tmTemp, Constants.EGM96_EARTH_MU),
				    		currentObj.prop.getAttitude(tmTemp, pv), odcfg.SpaceObject.Mass);


				    
				    if (combmeas)
				    {
						double[] fitv = currentObj.associatedObs.get(measNum).estimate(1, 1, sstaTemp).getEstimatedValue();
						spupd.setColumn(i, fitv);
						

				    }
				    else
				    {
						double[] fitv = currentObj.associatedObs.get(2 * measNum).estimate(1, 1, sstaTemp).getEstimatedValue();
						spupd.setEntry(0, i, fitv[0]);

					if (currentObj.Rsize > 1)
					{
					    fitv = fitv = currentObj.associatedObs.get(2 * measNum + 1).estimate(1, 1, sstaTemp).getEstimatedValue();
					    spupd.setEntry(1, i, fitv[0]);
					}

				    }
				    

				    
				}

				RealMatrix PyyTemp = currentObj.R.copy();
				RealVector yhatpreTemp = addColumns(spupd).mapMultiplyToSelf(currentObj.weight);
				for (int i = 0; i < numsig; i++)
				{
				    RealVector y = spupd.getColumnVector(i).subtract(yhatpreTemp);
				    PyyTemp = PyyTemp.add(y.outerProduct(y).scalarMultiply(currentObj.weight));
				}
				
				currentObj.objResults.Estimation.get(measNum).InnovationCovariance = PyyTemp.getData();;
				
				
				
			}

			

			//compute RIC covar
		    
		    for(int i = 0; i < currentObj.objResults.Estimation.size(); i ++)
		    {
		    	
		    	RealMatrix posTemp = new Array2DRowRealMatrix(Arrays.copyOfRange(currentObj.objResults.Estimation.get(i).EstimatedState,0,3));
		    	RealMatrix velTemp = new Array2DRowRealMatrix(Arrays.copyOfRange(currentObj.objResults.Estimation.get(i).EstimatedState,3,6));

		    	RealMatrix covarTemp = new Array2DRowRealMatrix(currentObj.objResults.Estimation.get(i).EstimatedCovariance);
		    	
		    	
		    	RealMatrix hhat = crossProduct(posTemp, velTemp);
		    	
		    	RealMatrix Rhat = posTemp.scalarMultiply(1 / posTemp.getFrobeniusNorm());
		    	RealMatrix Chat = (hhat).scalarMultiply(1 / (hhat).getFrobeniusNorm());
		    	RealMatrix Ihat = crossProduct(Chat,Rhat);
		    	

		    	
		    	
		    	RealMatrix RotationMat = new Array2DRowRealMatrix(3,3);
		    	
		    	RotationMat.setRowMatrix(0,Rhat.transpose());

		    	RotationMat.setRowMatrix(1,Ihat.transpose());

		    	RotationMat.setRowMatrix(2,Chat.transpose());	    	
		    	
		    	currentObj.objResults.Estimation.get(i).RICCovariance = RotationMat.multiply(covarTemp.getSubMatrix(0,2,0,2)).multiply(RotationMat.transpose()).getData();

		    }

			results = new JSONResults();
			
			results.Estimation.addAll(currentObj.objResults.Estimation);
		    results.Propagation.Time = currentObj.objResults.Propagation.Time;
		    results.Propagation.State = currentObj.objResults.Propagation.State;
		    


		    
		    resultsArr.add(results);
		}
	
		
}

	ArrayList<Hypothesis> GenerateHypotheses(Measurements.JSONMeasurement obs, GroundStation station, Frame propframe, AbsoluteDate date, double sigmaRA, double sigmaRAd, double sigmaDec, double sigmaDecd)
	{

		double RA = obs.RightAscension;
		double RA_d = obs.RightAscensionRate;
		
		double Dec = obs.Declination;
		double Dec_d = obs.DeclinationRate;
		
		
		TimeStampedPVCoordinates stationCoords = station.getBaseFrame().getPVCoordinates(date, propframe);
		
		
		ArrayList <CAR.CARGaussianElement> CARGaussians = new CAR(RA, RA_d, Dec, Dec_d, stationCoords, 20000.0, 1.0, 5000.0 ,17000000, 22000000.0, 0.04).getCAR();

		ArrayList<Hypothesis> objectHypotheses= new ArrayList<Hypothesis>();
		
		
		System.out.println(CARGaussians.size());

		
		for(int i = 0; i < CARGaussians.size(); i++)
		{
			System.out.println("LT Adjust:" + i);
			Hypothesis singleHypothesis= new Hypothesis();

			
			RealVector meanTemp = new ArrayRealVector(new double[] {CARGaussians.get(i).rangeMean, RA, Dec, CARGaussians.get(i).rangeRateMean, RA_d, Dec_d });
			RealMatrix CovarTemp = new DiagonalMatrix(new double[] {CARGaussians.get(i).rangeStd, sigmaRA, sigmaDec, CARGaussians.get(i).rangeRateStd, sigmaRAd, sigmaDecd });
			
			RealMatrix sigma = GenerateSigmaPoints(meanTemp, CovarTemp, 6, 6*2);
			RealMatrix sigmaXYZ = new Array2DRowRealMatrix(6,6*2);

			
			
			for (int j = 0; j < sigmaXYZ.getColumnDimension(); j++)
			{
				
				double[] Xout = RangeRaDec2XYZ(sigma.getColumnVector(j), date, station, propframe);
		

				//adjust for LT
				
			    Xout = new ManualPropagation(promotedTracks.get(0).odcfg, 6).propagate(0, Arrays.copyOfRange(Xout, 0, 6), sigma.getEntry(0,j)/Constants.SPEED_OF_LIGHT);
				
				sigmaXYZ.setColumn(j, Xout);

			}

			
			
			RealMatrix Pxx = new Array2DRowRealMatrix(6, 6);
			RealVector xhat = addColumns(sigmaXYZ).mapMultiplyToSelf(0.5/6);
			
			for (int j = 0; j < sigmaXYZ.getColumnDimension(); j++)
			{
			    RealVector xhatdiff = sigmaXYZ.getColumnVector(j).subtract(xhat);
			    Pxx = Pxx.add(xhatdiff.outerProduct(xhatdiff).scalarMultiply(0.5/6));
			}
			
			
			//return mean/covar structure
			singleHypothesis.xhat = new ArrayRealVector(xhat.toArray());
			singleHypothesis.P = Pxx;
			singleHypothesis.weight = CARGaussians.get(i).weight;

			
			objectHypotheses.add(singleHypothesis);
			
			
			
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		
		try {
			FileWriter writer = new FileWriter("../../../Matlab/OriginaloutputStates.txt"); 

		for(Hypothesis line: objectHypotheses) {
			
			for(int i = 0; i <  line.xhat.getDimension(); i++)
			{
				writer.write(line.xhat.getEntry(i)+",");
			}
			
			for(int i = 0; i <  line.P.getRowDimension(); i++)
			{
				for(int j = 0; j <  line.P.getColumnDimension(); j++)
				{
					writer.write(line.P.getEntry(i,j)+",");
				}
			}
			
			
		  writer.write(line.weight+"," + date + System.lineSeparator());
		  }
		
		writer.close();
	
		} catch (IOException e) {
	
		e.printStackTrace();
		System.exit(0);
		} 
		
		/////////////////////////////////////////////////////////////////////////////////////////
		
		
		return objectHypotheses;
	}
	
	
	
	
	RealMatrix GenerateSigmaPoints(RealVector xhat, RealMatrix P, int numsta, int numsig)
	{
		RealMatrix sigma = new Array2DRowRealMatrix(numsta, numsig);
		
		RealMatrix Ptemp = P.scalarMultiply(numsta);
		
		RealMatrix sqrP = new CholeskyDecomposition(
				Ptemp.add(Ptemp.transpose()).scalarMultiply(0.5), 1E-6, 1E-14).getL();

		for (int j = 0; j < numsta; j++)
		{
			sigma.setColumnVector(j, xhat.add(sqrP.getColumnVector(j)));
			sigma.setColumnVector(numsta + j, xhat.subtract(sqrP.getColumnVector(j)));
		}
		
		
		return sigma;
	}

	
	
	
	
	
	double[] RangeRaDec2XYZ(RealVector RangeRaDec, AbsoluteDate date, GroundStation stat, Frame inertialFrame)
	{
		
		//CAR FrameTransfer Test
		
		double Range = RangeRaDec.getEntry(0);
		double Range_d = RangeRaDec.getEntry(3);		
		double RA = RangeRaDec.getEntry(1);
		double RA_d = RangeRaDec.getEntry(4);		
		double Dec = RangeRaDec.getEntry(2);
		double Dec_d = RangeRaDec.getEntry(5);
		
		
		//Compute Inertial pos/vel relative to station
		Vector3D topoPos = new Vector3D(new double[] {Range*Math.cos(Dec)*Math.cos(RA), Range*Math.cos(Dec)*Math.sin(RA), Range*Math.sin(Dec)});
		Vector3D topoVel = new Vector3D(new double[] {Range_d*Math.cos(Dec)*Math.cos(RA) - Range*Math.sin(Dec)*Math.cos(RA)*Dec_d - Range*Math.cos(Dec)*Math.sin(RA)*RA_d,
															   Range_d*Math.cos(Dec)*Math.sin(RA) - Range*Math.sin(Dec)*Math.sin(RA)*Dec_d + Range*Math.cos(Dec)*Math.cos(RA)*RA_d,
															   Range_d*Math.sin(Dec) + Range*Math.cos(Dec) * Dec_d});
		

		//get station coords
		//add station coords to observation coords.
		Vector3D Pos = new Vector3D(1, topoPos, 1, stat.getBaseFrame().getPVCoordinates(date, inertialFrame).getPosition());
		Vector3D Vel = new Vector3D(1, topoVel, 1, stat.getBaseFrame().getPVCoordinates(date, inertialFrame).getVelocity());

		
		
		return new double[] {Pos.getX(), Pos.getY(), Pos.getZ(), Vel.getX(), Vel.getY(), Vel.getZ()};
				
		
	}
	
	
	
	RealVector quadCheck(RealVector measurement, RealVector state)
	{
	    //Innovations QuadCheck
		RealVector Innov = measurement.subtract(state);
	
		if(combmeas && Innov.getEntry(0) > Math.PI)
		{
			Innov.setEntry(0, Innov.getEntry(0) - 2 * Math.PI);
		}
		else if(combmeas && Innov.getEntry(0) < -1*  Math.PI)
		{
			Innov.setEntry(0, Innov.getEntry(0) + 2 * Math.PI);
		}
		
		return Innov;
	}
	
	
	
	
	
	RealMatrix crossProduct(RealMatrix vect_A, RealMatrix vect_B) 
	{ 
		RealMatrix cross_P = new Array2DRowRealMatrix(3,1);
		
	    cross_P.setEntry(0,0, vect_A.getEntry(1,0) * vect_B.getEntry(2,0) - vect_A.getEntry(2,0) * vect_B.getEntry(1,0)); 
	    cross_P.setEntry(1,0, vect_A.getEntry(2,0) * vect_B.getEntry(0,0) - vect_A.getEntry(0,0) * vect_B.getEntry(2,0)); 
	    cross_P.setEntry(2,0, vect_A.getEntry(0,0) * vect_B.getEntry(1,0) - vect_A.getEntry(1,0) * vect_B.getEntry(0,0)); 
	    
	    return cross_P;
	}
	
	
	ArrayList<ArrayList<JPDALikelihoods>> JPDAJointEvents(ArrayList<ArrayList<JPDALikelihoods>> JointEvents, double[][] MarginalEvents, ArrayList<JPDALikelihoods> SingleJointEvent, int objNum)
	{
		
		for(int measNum = 0; measNum<MarginalEvents[objNum].length; measNum++)
		{	
			ArrayList<JPDALikelihoods> SingleJointEventTemp = new ArrayList<JPDALikelihoods>();
			//Copy Single Joint Event into Single Joint Event Temp
			for(int m = 0; m < SingleJointEvent.size(); m++)
			{
				SingleJointEventTemp.add(SingleJointEvent.get(m));
			}
			//Add event if that measurement has not been used
			boolean skip = false;
			for(int m = 0; m < SingleJointEvent.size(); m++)
			{
				if(measNum == SingleJointEventTemp.get(m).measurement && measNum != 0)
				{
					skip = true;
				}
			}
	
			if(skip == true)
			{
				continue;
			}
			
			JPDALikelihoods temp = new JPDALikelihoods();
			
			temp.measurement = measNum;
			temp.object = objNum;
			temp.Psi = MarginalEvents[objNum][measNum];
			
			SingleJointEventTemp.add(temp);
			//decide to repeat or not
			if(MarginalEvents.length == objNum+1)
			{
				JointEvents.add(SingleJointEventTemp);
			}
			else
			{
				JointEvents = JPDAJointEvents(JointEvents, MarginalEvents, SingleJointEventTemp, objNum+1);
			}
		}
		

		return JointEvents;

	}
	
	double[] stack(RealMatrix mat, double[] arr)
	{
	    int i,j;
	    int m = mat.getRowDimension();
	    int n = mat.getColumnDimension();
	    double[][] matdata = mat.getData();

	    for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		    arr[i*m + j] = matdata[j][i];

	    return(arr);
	}

	RealMatrix unstack(RealMatrix mat, double[] arr)
	{
	    int m = mat.getRowDimension();
	    int n = mat.getColumnDimension();

	    for (int i = 0; i < n; i++)
		mat.setColumn(i, Arrays.copyOfRange(arr, i*m, (i+1)*m));

	    return(mat);
	}

	ArrayRealVector addColumns(RealMatrix mat)
	{
	    int i,j;
	    double sum;
	    int m = mat.getRowDimension();
	    int n = mat.getColumnDimension();
	    double[][] arr = mat.getData();
	    ArrayRealVector out = new ArrayRealVector(m);

	    for (j = 0; j < m; j++)
	    {
		sum = 0.0;
		for (i = 0; i < n; i++)
		    sum += arr[j][i];
		out.setEntry(j, sum);
	    }

	    return(out);
	}
    }

    
    class Hypothesis
    {
    	
        ArrayRealVector xhat;
 	    RealMatrix P;    
 	    double weight;
 	    
    }
    
    
    class JSONResults
    {
		class JSONEstimation
		{
		    String Time;
		    HashMap<String, double[]> PreFit;
		    HashMap<String, double[]> PostFit;
		    double[] EstimatedState;
		    double[] EstimatedAcceleration;
		    double[][] EstimatedCovariance;
		    double[][] InnovationCovariance;
		    double[][] RICCovariance;
		    int numConsideredParams;
	
		    
		    
		    public JSONEstimation()
		    {
			PreFit = new HashMap<String, double[]>();
			PostFit = new HashMap<String, double[]>();
			
		    }
		}
	
		class JSONPropagation
		{
		    String Time;
		    double[] State;
		}
	
		String Filter;
		ArrayList<JSONEstimation> Estimation;
		JSONPropagation Propagation;
	
		public JSONResults()
		{
		    Estimation = new ArrayList<JSONEstimation>();
		    Propagation = new JSONPropagation();
		}
    }
    
    

	class JPDALikelihoods
	{
		double Psi;
		int object;
		int measurement;
		
	}


    
    class SCstate
    {
        Settings odcfg;

	    SpacecraftState[] ssta;
	    ManualPropagation prop;
	    
        AbsoluteDate tm;
        AbsoluteDate epoch;

	    double[] Xi;
	    double[][] Pi;
	    RealVector xhatpre;
        ArrayRealVector xhat;
	    RealMatrix Ppre;
	    RealMatrix P;
		Array2DRowRealMatrix sigpr;
	    RealMatrix Q;
	    Array2DRowRealMatrix R;
		RealVector yhatpre;
	    RealMatrix Pyy;
		RealMatrix K;
		
		double hypothesisWeight;
		
	    int Rsize;
	    double ProbD;
	    double weight;
	    double[] spvec;
	    
    	ArrayList<JPDALikelihoods> MarginalEvents;

	    
		boolean McReynoldsConsistencyPass;
		boolean dataAssociated;
		boolean badHypothesis;
		
		smootherData smootherResults;
		JSONResults objResults;
		JSONResults.JSONEstimation odout;
		ArrayList<Measurements.JSONMeasurement> associatedObsJSON;
		ArrayList<ObservedMeasurement> associatedObs;

		Double RaNoise;
		Double DecNoise;
		
	    class smootherData
	    {
	    	class smootherStep
	    	{
	    		
	    		RealMatrix Ppre;
	    		RealMatrix Ppost;
	    		RealMatrix xpre;
	    		RealMatrix xpost;
	    		RealMatrix sigPre;
	    		RealMatrix sigPost;
	    		
	    		RealMatrix xstar;
	    		RealMatrix Pstar;
	    		
	    		AbsoluteDate tmSmoother;
	    		ObservedMeasurement measObjsSmoother;
	    		ObservedMeasurement measObjsSmootherNoComb;

	    	}
	    	
	    	ArrayList<smootherStep> smoother;
	    	
	    	public smootherData()
	    	{
	    		smoother = new ArrayList<smootherStep>();
	    	}
	    }	
	    
	    
	    public SCstate()
	    {
			objResults = new JSONResults();
			associatedObsJSON = new ArrayList<Measurements.JSONMeasurement>();
			associatedObs = new ArrayList<ObservedMeasurement>();
			MarginalEvents = new ArrayList<JPDALikelihoods>();


	    }
		
    }
    
}
