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

public class Estimation
{
    protected Settings odcfg;
    protected Measurements odobs;

    protected String[] meanames;
    protected boolean combmeas;

    protected JSONResults results;
    protected JSONResults[] resultsArr;

    

    public final static String DMC_ACC_ESTM = "DMCaccest";
    public final static String DMC_ACC_PROP = "DMCaccprop";

    
	ArrayList<SCstate> states = new ArrayList<SCstate>();
    
    
    public Estimation(String[] cfgjson, String obsjson)
    {

    for(int i = 0; i < cfgjson.length; i++)
    {

		states.add(i, new SCstate());

		states.get(i).odcfg = Settings.loadJSON(cfgjson[i]);

		if (states.get(i).odcfg.Estimation.Filter == null)
			states.get(i).odcfg.Estimation.Filter = "UKF";

		if (states.get(i).odcfg.Estimation.Filter.equals("UKF") && states.get(i).odcfg.Gravity.Degree >= 2 && states.get(i).odcfg.Gravity.Order >= 0)
			states.get(i).odcfg.forces.add(0, new NewtonianAttraction(Constants.EGM96_EARTH_MU));

		if (states.get(i).odcfg.Estimation.NoiseTimeDelta <= 0.0)
			states.get(i).odcfg.Estimation.NoiseTimeDelta = 10.0;

    }
		odobs = Measurements.loadJSON(states.get(0).odcfg, obsjson);

		meanames = states.get(0).odcfg.Measurements.keySet().toArray(new String[0]);
		combmeas = meanames[0].equals("Azimuth") || meanames[0].equals("Elevation") ||
		    meanames[0].equals("RightAscension") || meanames[0].equals("Declination") ||
		    meanames[0].equals("PositionVelocity");

    }

    public String[] determineOrbit()
    {

	results = new JSONResults();
	results.Filter = states.get(0).odcfg.Estimation.Filter;

	if (states.get(0).odcfg.Estimation.Filter.equals("UKF"))
	    new UnscentedKalmanFilter().determineOrbit();
	else
	    new ExtendedKalmanFilter().determineOrbit();


	String[] JSONOutput = new String[2*resultsArr.length];

	for(int i = 0; i < resultsArr.length; i++)
	{
		JSONOutput[i] =new GsonBuilder().setPrettyPrinting().create().toJson(resultsArr[i]);
		JSONOutput[i+resultsArr.length] =new GsonBuilder().setPrettyPrinting().create().toJson(states.get(i).associatedObs);

	}

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
	    int numsta = states.get(0).odcfg.estparams.size() + 6;
	    int numsig = 2*numsta;
	    int veclen = numsta*numsig;
		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			Settings odcfg = states.get(objNum).odcfg;
			states.get(objNum).P = new DiagonalMatrix(odcfg.Estimation.Covariance);
			states.get(objNum).Q = odcfg.getProcessNoiseMatrix();
			
		    states.get(objNum).Rsize = 0;
		    for (String s: meanames)
	    	states.get(objNum).Rsize += odcfg.Measurements.get(s).Error.length;

		    states.get(objNum).R = new Array2DRowRealMatrix(states.get(objNum).Rsize, states.get(objNum).Rsize);
		    for (int i = 0, j = 0; i < meanames.length; i++)
		    {
			Settings.JSONMeasurement jm = odcfg.Measurements.get(meanames[i]);
			for (int k = 0; k < jm.Error.length; k++)
			{
				states.get(objNum).R.setEntry(j, j, jm.Error[k]*jm.Error[k]);
			    j++;
			}
		    }

		    states.get(objNum).Xi = odcfg.getInitialState();
		    states.get(objNum).epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.Start),
							  DataManager.utcscale);
	
		    if (odobs.rawmeas.length > 0)
		    {
			AbsoluteDate tmto = new AbsoluteDate(DateTimeComponents.parseDateTime(odobs.rawmeas[0].Time),
							     DataManager.utcscale);
			if (Math.abs(tmto.durationFrom(states.get(objNum).epoch)) > 1.0)
			{
			    ManualPropagation prop0 = new ManualPropagation(odcfg, 6);
			    double[] Xout = prop0.propagate(0, Arrays.copyOfRange(states.get(objNum).Xi, 0, 6), tmto.durationFrom(states.get(objNum).epoch));
			    states.get(objNum).epoch = new AbsoluteDate(tmto, 0.0);
			    System.arraycopy(Xout, 0, states.get(objNum).Xi, 0, 6);
			}
		    }
	


		    states.get(objNum).xhat = new ArrayRealVector(states.get(objNum).Xi);
		    states.get(objNum).xhatpre = new ArrayRealVector(states.get(objNum).Xi);
		    states.get(objNum).weight = 0.5/numsta;
		    states.get(objNum).spvec = new double[veclen];

		    states.get(objNum).tm = new AbsoluteDate(states.get(objNum).epoch, 0.0);
		    states.get(objNum).ssta = new SpacecraftState[1];
	
		    states.get(objNum).prop = new ManualPropagation(odcfg, veclen);
		    
		    states.get(objNum).ProbD = odcfg.Estimation.ProbDetection;
			states.get(objNum).smootherResults = states.get(objNum).new smootherData();

		}
		
		for(int smootherIter = 0; smootherIter < states.get(0).odcfg.Estimation.USKFIterations; smootherIter++)
		{
			
			if(smootherIter >= 1)
			{
				for(int objNum = 0; objNum < states.size(); objNum++)
				{
					states.get(objNum).xhat = new ArrayRealVector(states.get(objNum).objResults.Estimation.get(0).EstimatedState);
			    	states.get(objNum).xhatpre = new ArrayRealVector(states.get(objNum).objResults.Estimation.get(0).EstimatedState);
			    	states.get(objNum).P = new DiagonalMatrix(states.get(objNum).odcfg.Estimation.Covariance);
				    states.get(objNum).tm = new AbsoluteDate(states.get(objNum).objResults.Estimation.get(0).Time, DataManager.utcscale);
					states.get(objNum).smootherResults = states.get(objNum).new smootherData();
					
					states.get(objNum).objResults = new JSONResults();
					states.get(objNum).associatedObs = new ArrayList<Measurements.JSONMeasurement>();
					
					
					states.get(objNum).obsConsidered = 0;
				}
			}
			
    	System.out.printf("Iteration: %d \n", smootherIter);

			
		int mix = -1;
		int additionalMeas = 0;
		
		boolean finished = false;

		//// Loop for all obs starts
		EstimationLoop:
	    while(!finished)
	    {

	    mix = mix + additionalMeas + 1;
    	System.out.printf("Percent Complete: %.2f%% \n", (double) mix/odobs.rawmeas.length * 100);

	    //calc additional meas
	    additionalMeas = 0;

		if (mix + additionalMeas < odobs.rawmeas.length)
		{
		    while(true)
		    {	    	

		    	int i = 1; 
		    	if(mix+additionalMeas+1 < odobs.rawmeas.length && odobs.rawmeas[mix+additionalMeas].Time.equals(odobs.rawmeas[mix+additionalMeas+1].Time)  
		    			&& odobs.rawmeas[mix+additionalMeas].Station.equals(odobs.rawmeas[mix+additionalMeas+1].Station))
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

			finished = true;
		}
		

		
		JPDAEvent[] MarginalEvents = new JPDAEvent[states.size()];
		ArrayList<JPDAEvent> JointEvents = new ArrayList<JPDAEvent>();
		JPDAEvent SingleJointEvent = new JPDAEvent();
		

		///////////////////////////////////////////////////////// FIRST OBJECT LOOP /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////// Propagation & Constructing Marginal Events/////////////////////////////////

		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			
		SCstate currentObj = states.get(objNum);
		odcfg = currentObj.odcfg;
		
		currentObj.odout = results.new JSONEstimation();

		Array2DRowRealMatrix sigma = new Array2DRowRealMatrix(numsta, numsig);
		Array2DRowRealMatrix spupd = new Array2DRowRealMatrix(states.get(0).Rsize, numsig);
	    
	    
		AbsoluteDate t0 = currentObj.tm;
		
		if(!finished)
		{
			
			currentObj.tm = new AbsoluteDate(DateTimeComponents.parseDateTime(odobs.rawmeas[mix].Time),
					  DataManager.utcscale);

		    
		    double[] pv = currentObj.xhat.toArray();
		    TimeStampedPVCoordinates pvs = new TimeStampedPVCoordinates(currentObj.tm,
										new Vector3D(pv[0], pv[1], pv[2]),
										new Vector3D(pv[3], pv[4], pv[5]));

		    if (odcfg.stations != null && odobs.rawmeas[mix].Station != null)
		    {

			PVCoordinates pvi = odcfg.stations.get(odobs.rawmeas[mix].Station).getBaseFrame().
			    getPVCoordinates(currentObj.tm, odcfg.propframe);
			
		    }

		    states.get(objNum).odout.Time = odobs.rawmeas[mix].Time;

		}
		else
		{

			currentObj.tm = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.End),
					  DataManager.utcscale);
		    

		}

		
		if(states.get(objNum).dataAssociated == true || mix == 0)
		{
			
		states.get(objNum).dataAssociated = false;
		
		currentObj.sigpr = new Array2DRowRealMatrix(numsta, numsig);

		
		RealMatrix Ptemp = currentObj.P.scalarMultiply(numsta);
		RealMatrix sqrP = new CholeskyDecomposition(
		    Ptemp.add(Ptemp.transpose()).scalarMultiply(0.5), 1E-6, 1E-14).getL();

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
			sigma = states.get(objNum).sigpr;
			
			
		}

		
		double propt0 = t0.durationFrom(currentObj.epoch);
		double propt1 = currentObj.tm.durationFrom(currentObj.epoch);
		if (propt0 == propt1)
			currentObj.sigpr.setSubMatrix(sigma.getData(), 0, 0);
		else
		    unstack(currentObj.sigpr, currentObj.prop.propagate(propt0, stack(sigma, currentObj.spvec), propt1));

		
		currentObj.xhatpre = addColumns(currentObj.sigpr).mapMultiplyToSelf(currentObj.weight);


		
		if (finished) 
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
			if(MarginalEvents[objNum] == null)
			{
			    MarginalEvents[objNum] = new JPDAEvent();
				JPDAEvent.JPDALikelihoods JPDAtemp2 = MarginalEvents[objNum].new JPDALikelihoods();
				JPDAtemp2.Psi = 1 - currentObj.ProbD;
				JPDAtemp2.object = objNum;
				JPDAtemp2.measurement = 0;

				MarginalEvents[objNum].Likelihoods.add(JPDAtemp2);

				
			}
		


		// run through jpda algorithms and calculate psi_r_i
			
	    //Innovations QuadCheck
		RealVector Innov = raw.subtract(currentObj.yhatpre);
	
		if(combmeas && Innov.getEntry(0) > Math.PI)
		{
			Innov.setEntry(0, Innov.getEntry(0) - 2 * Math.PI);
		}
		else if(combmeas && Innov.getEntry(0) < -1*  Math.PI)
		{
			Innov.setEntry(0, Innov.getEntry(0) + 2 * Math.PI);
		}
							
		RealMatrix MahalaTemp = MatrixUtils.createColumnRealMatrix(Innov.toArray());

		RealMatrix JPDATemp = MahalaTemp.transpose().multiply(MatrixUtils.inverse(currentObj.Pyy).multiply(MahalaTemp));

		
		if( odcfg.Estimation.GatingThreshold > Math.sqrt(JPDATemp.getEntry(0,0)))
		{
	
			JPDAEvent.JPDALikelihoods JPDAtemp2 = MarginalEvents[objNum].new JPDALikelihoods();


			//JPDAtemp2.Psi =  Math.exp(-JPDATemp.getEntry(0,0)/2) / (Math.sqrt(new LUDecomposition(currentObj.Pyy).getDeterminant() * Math.pow(2 * Math.PI, currentObj.Rsize)));
			JPDAtemp2.Psi = 1 - new ChiSquaredDistribution(currentObj.Rsize).cumulativeProbability(JPDATemp.getEntry(0,0));
			
			
			
			JPDAtemp2.object = objNum;
			JPDAtemp2.measurement = measNum+1;
			
			MarginalEvents[objNum].Likelihoods.add(JPDAtemp2);
			

		}

					
		// end both loops (measurement then object)
		}


    		if(MarginalEvents[objNum].Likelihoods.size() > 1)
    		{
    			states.get(objNum).obsConsidered++;
    		}

        
	    }
		
	

		
		//Determine joint events
		
		JointEvents = JPDAJointEvents(JointEvents, MarginalEvents, SingleJointEvent);

		//compute Probabilities
		double[] JPDAProbability = new double[JointEvents.size()];
		Arrays.fill(JPDAProbability, 1);
		
		for(int i = 0; i<JointEvents.size(); i++)
		{
			for(int j = 0; j<JointEvents.get(i).Likelihoods.size(); j++)
			{
				//skip cases where  only one marginal event to choose from. This is important for when Pd = 1, and an object does not 
				//have an associated observation.
				
				//This fix should not affect non Pd = 1 cases, since there is only one event to choose from in the first place, and the 
				//probabilities get normalized after this.
				if(MarginalEvents[j].Likelihoods.size() != 1)
				{
					JPDAProbability[i] = JPDAProbability[i] * JointEvents.get(i).Likelihoods.get(j).Psi;
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


		for(int i = 0; i < JointEvents.get(maxProbIndex).Likelihoods.size(); i++)
		{
			if(JointEvents.get(maxProbIndex).Likelihoods.get(i).measurement != 0)
			{
				states.get(JointEvents.get(maxProbIndex).Likelihoods.get(i).object).associatedObs.add( 
						odobs.rawmeas[mix + JointEvents.get(maxProbIndex).Likelihoods.get(i).measurement - 1]);
				states.get(JointEvents.get(maxProbIndex).Likelihoods.get(i).object).dataAssociated = true;
			}
		
		}



		
		///////////////////////////////////////////////////////// SECOND OBJECT LOOP /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////// Compute probability of joint events, Perform update step ///////////////////

		for(int objNum = 0; objNum < states.size(); objNum++)
		{		
			SCstate currentObj = states.get(objNum);
			odcfg = currentObj.odcfg;
			
			double Beta = 0;
			RealMatrix Betav = new Array2DRowRealMatrix(currentObj.Rsize, 1);
			RealMatrix Betavvt = new Array2DRowRealMatrix(currentObj.Rsize, currentObj.Rsize);



			for(int measNum = 0; measNum <=additionalMeas; measNum++)
			{
			
				//compute Beta
				//check each probability to see if it has the corresponding obj/meas combo

				for(int i = 0; i < JointEvents.size(); i++)
				{

					for(int j = 0; j < JointEvents.get(i).Likelihoods.size(); j++)
					{

						if(JointEvents.get(i).Likelihoods.get(j).object == objNum && JointEvents.get(i).Likelihoods.get(j).measurement - 1 == measNum)
						{
							
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
						    RealMatrix Innovation = MatrixUtils.createColumnRealMatrix(raw.subtract(currentObj.yhatpre).toArray());

						    
				    		if(combmeas && Innovation.getEntry(0,0) > Math.PI)
							{
				    			Innovation.setEntry(0,0, Innovation.getEntry(0,0) - 2 * Math.PI);
							}
							else if(combmeas && Innovation.getEntry(0,0) < -1 *  Math.PI)
							{
								Innovation.setEntry(0,0, Innovation.getEntry(0,0) + 2 * Math.PI);
							}
							

							
							
							Beta = Beta + JPDAProbability[i];
							Betav= Betav.add(Innovation.scalarMultiply(JPDAProbability[i]));
							Betavvt = Betavvt.add(Innovation.multiply(Innovation.transpose()).scalarMultiply(JPDAProbability[i]));


						}
						
					}
					
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
					states.get(objNum).odout.PreFit.put(meanames[i], currentObj.yhatpre.toArray());
					states.get(objNum).odout.PostFit.put(meanames[i], fitv);
				}
				else
				{
					states.get(objNum).odout.PreFit.put(meanames[i], new double[] {currentObj.yhatpre.getEntry(i)});
					states.get(objNum).odout.PostFit.put(meanames[i], new double[] {fitv[i]});
				    
				    
				}
			    }
			}
			else
			{

			    double[] fitv = odobs.measobjs.get(mix*2).estimate(1, 1, currentObj.ssta).getEstimatedValue();
			    states.get(objNum).odout.PreFit.put(meanames[0], new double[] {currentObj.yhatpre.getEntry(0)});
			    states.get(objNum).odout.PostFit.put(meanames[0], fitv);
	
			    if (currentObj.Rsize > 1)
			    {
				fitv = odobs.measobjs.get(mix*2 + 1).estimate(1, 1, currentObj.ssta).getEstimatedValue();
				states.get(objNum).odout.PreFit.put(meanames[1], new double[] {currentObj.yhatpre.getEntry(1)});
				states.get(objNum).odout.PostFit.put(meanames[1], fitv);
			    }
			}

			if(states.get(objNum).dataAssociated == true)
			states.get(objNum).objResults.Estimation.add(states.get(objNum).odout);

			
			
			
			if(odcfg.Estimation.Smoother.equals("On") && states.get(objNum).dataAssociated == true)
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

				SCstate.smootherData.smootherStep smout = states.get(objNum).smootherResults.new smootherStep();
		
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
				
		    

		}
		

	    }
		
		///////////////////////////////////////////////////////// THIRD OBJECT LOOP /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////// Smooth Data ///////////////////////////////////////////////////////////////

		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			SCstate currentObj = states.get(objNum);
			odcfg = currentObj.odcfg;
			
		    double[] pv = currentObj.xhatpre.toArray();


		    currentObj.tm = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.End),
					  DataManager.utcscale);

		    CartesianOrbit cart  = new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										new Vector3D(pv[3], pv[4], pv[5])),
							      odcfg.propframe, currentObj.tm, Constants.EGM96_EARTH_MU);
		    if (cart.getA() <= Constants.WGS84_EARTH_EQUATORIAL_RADIUS)
			throw(new RuntimeException(String.format("Invalid semi-major axis %f", cart.getA())));

		    
		    
		    states.get(objNum).objResults.Propagation.Time = odcfg.Propagation.End;
		    states.get(objNum).objResults.Propagation.State = pv;
		    

		    
		    // run smoother
	
			if(odcfg.Estimation.Smoother.equals("On"))
			{
				
				states.get(objNum).McReynoldsConsistencyPass = true;
				
				
				int smSize = states.get(objNum).smootherResults.smoother.size()-1;

				states.get(objNum).smootherResults.smoother.get(smSize).xstar 
					= states.get(objNum).smootherResults.smoother.get(smSize).xpost;
				states.get(objNum).smootherResults.smoother.get(smSize).Pstar 
					= states.get(objNum).smootherResults.smoother.get(smSize).Ppost;			
			
			double McReynoldsVal = -99999;

		    for(int i = 0; i < states.get(objNum).smootherResults.smoother.size()-1;i++)
		    {

		    	SCstate.smootherData.smootherStep smDatak1 = states.get(objNum).smootherResults.smoother.get(smSize - i);
		    	SCstate.smootherData.smootherStep smDatak = states.get(objNum).smootherResults.smoother.get(smSize - i - 1);
	
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
				


				states.get(objNum).smootherResults.smoother.set(smSize - i - 1,smDatak);
				
				if(states.get(objNum).McReynoldsConsistencyPass == true)
				{
					
					RealMatrix delx = smDatak.xpost.subtract(smDatak.xstar);
							
					RealMatrix delP = smDatak.Ppost.subtract(smDatak.Pstar);
					
					for(int j = 0; j < 5; j++)
					{			

						McReynoldsVal = Math.max(McReynoldsVal, Math.abs(delx.getEntry(j,0)) / Math.sqrt(Math.abs(delP.getEntry(j,j))));
						

						
						if(Math.abs(delx.getEntry(j,0)) / Math.sqrt(Math.abs(delP.getEntry(j,j))) >= 3)
						{
							states.get(objNum).McReynoldsConsistencyPass = false;

						}
					}
				}
		    }
		    
		    		    
		    System.out.println("McReynolds Pass: " + states.get(objNum).McReynoldsConsistencyPass +  "   McReynoldsVal: " + McReynoldsVal);
    		System.out.println("Object " + objNum + " Observations Considered: " + states.get(objNum).obsConsidered);
    		System.out.println("Object " + objNum + " Observations Associated: " + states.get(objNum).associatedObs.size());


		    //store in results
		    
			    SpacecraftState[] smSsta = new SpacecraftState[1];		
			    
			    for(int j = 0; j < states.get(objNum).smootherResults.smoother.size();j++)
			    {
			    	
				    pv = states.get(objNum).smootherResults.smoother.get(j).xstar.getColumn(0);
				    
				    currentObj.tm = states.get(objNum).smootherResults.smoother.get(j).tmSmoother;
				    
				    smSsta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
												       new Vector3D(pv[3], pv[4], pv[5])), odcfg.propframe, currentObj.tm, Constants.EGM96_EARTH_MU));
		
				    
			    	//compute smoothed residuals
				    
				    if (combmeas)
					{
					    for (int i = 0; i < meanames.length; i++)
					    {
						double[] fitv = states.get(objNum).smootherResults.smoother.get(j).measObjsSmoother.estimate(1, 1, smSsta).getEstimatedValue();
						if (meanames.length == 1)
						{

							states.get(objNum).objResults.Estimation.get(j).PostFit.put(meanames[i], fitv);
						}
						else
						{

							states.get(objNum).objResults.Estimation.get(j).PostFit.put(meanames[i], new double[] {fitv[i]});
							
		
						}
					    }
					}
					else
					{
					    double[] fitv = states.get(objNum).smootherResults.smoother.get(j).measObjsSmoother.estimate(1, 1, smSsta).getEstimatedValue();
					    states.get(objNum).objResults.Estimation.get(j).PostFit.put(meanames[0], fitv);
		
					    if (currentObj.Rsize > 1)
					    {

						fitv = states.get(objNum).smootherResults.smoother.get(j).measObjsSmootherNoComb.estimate(1, 1, smSsta).getEstimatedValue();
						states.get(objNum).objResults.Estimation.get(j).PostFit.put(meanames[1], fitv);
					    }
					}
			    	
					
					
					
			    	//store
					states.get(objNum).objResults.Estimation.get(j).EstimatedState = states.get(objNum).smootherResults.smoother.get(j).xstar.getColumn(0);
					states.get(objNum).objResults.Estimation.get(j).EstimatedCovariance = states.get(objNum).smootherResults.smoother.get(j).Pstar.getData();
			    		    	
			    }
			    
			}
		
		}

			//break out of smoother iteration if all data has been established
			boolean allDataAssociated = true;
			for(int objNum = 0; objNum < states.size(); objNum++)
			{
				if(states.get(objNum).McReynoldsConsistencyPass == false || states.get(objNum).obsConsidered != states.get(objNum).associatedObs.size())
				{
					allDataAssociated = false;
				}
						
			}
			
			if(allDataAssociated == true)
			{
				break;
			}
		
		}
		
		resultsArr = new JSONResults[states.size()];
		

		for(int objNum = 0; objNum < states.size(); objNum++)
		{
			
		    //compute RIC covar
		    
		    for(int i = 0; i < states.get(objNum).objResults.Estimation.size(); i ++)
		    {
		    	
		    	RealMatrix posTemp = new Array2DRowRealMatrix(Arrays.copyOfRange(states.get(objNum).objResults.Estimation.get(i).EstimatedState,0,3));
		    	RealMatrix velTemp = new Array2DRowRealMatrix(Arrays.copyOfRange(states.get(objNum).objResults.Estimation.get(i).EstimatedState,3,6));

		    	RealMatrix covarTemp = new Array2DRowRealMatrix(states.get(objNum).objResults.Estimation.get(i).EstimatedCovariance);
		    	
		    	
		    	RealMatrix hhat = crossProduct(posTemp, velTemp);
		    	
		    	RealMatrix Rhat = posTemp.scalarMultiply(1 / posTemp.getFrobeniusNorm());
		    	RealMatrix Chat = (hhat).scalarMultiply(1 / (hhat).getFrobeniusNorm());
		    	RealMatrix Ihat = crossProduct(Chat,Rhat);
		    	

		    	
		    	
		    	RealMatrix RotationMat = new Array2DRowRealMatrix(3,3);
		    	
		    	RotationMat.setRowMatrix(0,Rhat.transpose());

		    	RotationMat.setRowMatrix(1,Ihat.transpose());

		    	RotationMat.setRowMatrix(2,Chat.transpose());	    	
		    	
		    	states.get(objNum).objResults.Estimation.get(i).RICCovariance = RotationMat.multiply(covarTemp.getSubMatrix(0,2,0,2)).multiply(RotationMat.transpose()).getData();

		    }
			
			results = new JSONResults();
			
			results.Estimation.addAll(states.get(objNum).objResults.Estimation);
		    results.Propagation.Time = states.get(objNum).objResults.Propagation.Time;
		    results.Propagation.State = states.get(objNum).objResults.Propagation.State;
		    

		    
		    
		    
		    
		    resultsArr[objNum] = results;
		}

}

	
	RealMatrix crossProduct(RealMatrix vect_A, RealMatrix vect_B) 
	{ 
		RealMatrix cross_P = new Array2DRowRealMatrix(3,1);
		
	    cross_P.setEntry(0,0, vect_A.getEntry(1,0) * vect_B.getEntry(2,0) - vect_A.getEntry(2,0) * vect_B.getEntry(1,0)); 
	    cross_P.setEntry(1,0, vect_A.getEntry(2,0) * vect_B.getEntry(0,0) - vect_A.getEntry(0,0) * vect_B.getEntry(2,0)); 
	    cross_P.setEntry(2,0, vect_A.getEntry(0,0) * vect_B.getEntry(1,0) - vect_A.getEntry(1,0) * vect_B.getEntry(0,0)); 
	    
	    return cross_P;
	}
	
	
	ArrayList<JPDAEvent> JPDAJointEvents(ArrayList<JPDAEvent> JointEvents, JPDAEvent[] MarginalEvents, JPDAEvent SingleJointEvent)
	{
		for(int i = 0; i<MarginalEvents[0].Likelihoods.size(); i++)
		{	
			JPDAEvent SingleJointEventTemp = new JPDAEvent();
			//Copy Single Joint Event into Single Joint Event Temp
			for(int m = 0; m < SingleJointEvent.Likelihoods.size(); m++)
			{
				SingleJointEventTemp.Likelihoods.add(SingleJointEvent.Likelihoods.get(m));
			}
			//Add event if that measurement has not been used
			boolean skip = false;
			for(int m = 0; m < SingleJointEvent.Likelihoods.size(); m++)
			{
				if(MarginalEvents[0].Likelihoods.get(i).measurement == SingleJointEventTemp.Likelihoods.get(m).measurement && MarginalEvents[0].Likelihoods.get(i).measurement != 0)
				{
					skip = true;
				}
			}
	
			if(skip == true)
			{
				continue;
			}
			
			SingleJointEventTemp.Likelihoods.add(MarginalEvents[0].Likelihoods.get(i));
			//decide to repeat or not
			if(MarginalEvents.length == 1)
			{
				JointEvents.add(SingleJointEventTemp);
			}
			else
			{
				JointEvents = JPDAJointEvents(JointEvents, Arrays.copyOfRange(MarginalEvents,1,MarginalEvents.length), SingleJointEventTemp);
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
    
    
    class JPDAEvent
    {
    	
    	ArrayList<JPDALikelihoods> Likelihoods;
    	
    	class JPDALikelihoods
    	{
    		double Psi;
    		int object;
    		int measurement;  		
    		
    	}

    	public JPDAEvent()
    	{
    		Likelihoods = new ArrayList<JPDALikelihoods>();
    	}
    	
    }
    
    class SCstate
    {
    	int obsConsidered = 0; 
        Settings odcfg;

	    SpacecraftState[] ssta;
	    ManualPropagation prop;
	    
        AbsoluteDate tm;
        AbsoluteDate epoch;

	    double[] Xi;
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
		
	    int Rsize;
	    double ProbD;
	    double weight;
	    double[] spvec;
	    
		boolean McReynoldsConsistencyPass;
		boolean dataAssociated;

		smootherData smootherResults;
		JSONResults objResults;
		JSONResults.JSONEstimation odout;
		ArrayList<Measurements.JSONMeasurement> associatedObs;

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
			associatedObs = new ArrayList<Measurements.JSONMeasurement>();


	    }
		
    }
    
}
