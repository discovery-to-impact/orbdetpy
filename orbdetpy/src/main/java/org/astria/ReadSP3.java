/*
 * Estimation.java - Implementation of estimation algorithms.
 * Copyright (C) 2018-2020 University of Texas
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


import java.lang.Math;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.FileWriter;
import java.nio.file.Files;
import java.io.File;
import java.nio.file.StandardOpenOption;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.files.sp3.SP3Parser;
import org.orekit.files.sp3.SP3File;
import org.orekit.files.sp3.SP3File.SP3Ephemeris;
import org.orekit.files.ccsds.TDMFile;
import org.orekit.files.ccsds.TDMParser;
import org.orekit.files.ccsds.TDMFile.ObservationsBlock;
import org.orekit.files.ccsds.TDMParser.TDMFileFormat;
import org.orekit.propagation.BoundedPropagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.time.DateComponents;
import org.orekit.time.DateTimeComponents;
import org.orekit.time.GNSSDate;
import org.orekit.gnss.SatelliteSystem;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.frames.TopocentricFrame;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Transform;
import org.orekit.frames.ITRFVersion;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.models.earth.EarthStandardAtmosphereRefraction;
import org.orekit.estimation.measurements.modifiers.AngularRadioRefractionModifier;
import org.orekit.files.ccsds.CCSDSFrame;

public final class ReadSP3
{

    public ReadSP3()
    {

    	final ObservableSatellite obssat = new ObservableSatellite(0);
    	final double[] zeros = new double[]{0.0, 0.0};
    	final double[] ones = new double[]{1.0, 1.0};
    	
    	// NORAD ID to PRN data.
        HashMap<String, String> NORAD2PRN = new HashMap<String, String>();
        
        NORAD2PRN.put("2011-036A","G01");
        NORAD2PRN.put("2004-045A","G02");
        NORAD2PRN.put("2014-068A","G03");
        NORAD2PRN.put("2018-109A","G04");
        NORAD2PRN.put("2009-043A","G05");
        NORAD2PRN.put("2014-026A","G06");
        NORAD2PRN.put("2008-012A","G07");
        NORAD2PRN.put("2015-033A","G08");
        NORAD2PRN.put("2014-045A","G09");
        NORAD2PRN.put("2015-062A","G10");
        NORAD2PRN.put("1999-055A","G11");
        NORAD2PRN.put("2006-052A","G12");
        NORAD2PRN.put("1997-035A","G13");
        NORAD2PRN.put("2000-071A","G14");
        NORAD2PRN.put("2007-047A","G15");
        NORAD2PRN.put("2003-005A","G16");
        NORAD2PRN.put("2005-038A","G17");
        NORAD2PRN.put("1993-068A","G18");
        NORAD2PRN.put("2004-009A","G19");
        NORAD2PRN.put("2000-025A","G20");
        NORAD2PRN.put("2003-010A","G21");
        NORAD2PRN.put("2003-058A","G22");
        NORAD2PRN.put("2004-023A","G23");
        NORAD2PRN.put("2012-053A","G24");
        NORAD2PRN.put("2010-022A","G25");
        NORAD2PRN.put("2015-013A","G26");
        NORAD2PRN.put("2013-023A","G27");
        NORAD2PRN.put("2000-040A","G28");
        NORAD2PRN.put("2007-062A","G29");
        NORAD2PRN.put("2014-008A","G30");
        NORAD2PRN.put("2006-042A","G31");
        NORAD2PRN.put("2016-007A","G32");
        
        
        //Station lookup, station coords in degrees
        HashMap<String, double[]> StationLookup = new HashMap<String, double[]>();

        //Comments show groundstation coordinate reference, then reference frame for measurements.
        
        StationLookup.put("ART",new double[] {38.21586, -6.62771, 570}); //WGS, GCRF
        StationLookup.put("SMART-01-A-SUTH",new double[] {-32.38072, 20.81078, 1761}); //ITRF, EME2000
        StationLookup.put("SMART-01-B-SUTH",new double[] {-32.38072, 20.81078, 1761}); //ITRF, EME2000
        StationLookup.put("TRCK1",new double[] {38.54347, -4.40849, 1120}); //WGS, EME2000
        StationLookup.put("TRCK2",new double[] {38.54347, -4.40841, 1120}); //WGS, EME2000
        StationLookup.put("AROAC-T30",new double[] {46.71043251, 23.59323554, 783.405}); //WGS, ICRF
        StationLookup.put("BITNET-T30",new double[] {46.6760048, 23.1189174, 1170}); //WGS, ICRF
        StationLookup.put("MAM",new double[] {48.213256, 11.1766916, 582}); //WGS, ICRF
        StationLookup.put("PIRATE",new double[] {28.2992333, -16.51020277, 2422}); //WGS, ICRF
        StationLookup.put("COAST",new double[] {28.2990611, -16.5101778, 2420}); //WGS, ICRF
        
        AbsoluteDate referenceDate = null;
        System.out.println("adsasd");
        // read in TDM file using java.utils.
        File dir = new File("D:\\School\\Work\\NATOData\\TDM");
        File[] directoryListing = dir.listFiles();
        Arrays.sort(directoryListing);

		for (File child : directoryListing) {
			
		System.out.println(child);	
			
    	String TDMFile = child.toString();

    	ArrayList<ArrayList<Measurements.SimulatedMeasurement>> TDMData = Utilities.importTDM(TDMFile, "KEYVALUE");
    			
    	
    	// Extract Norad ID and convert to PRN
    	//System.out.println(TDMData.get(0).get(0).station);
    	String PRN = NORAD2PRN.get(TDMData.get(0).get(0).NORAD);
    	
    	if(PRN == null)
    	{
    		System.out.println("Unknown Participant");
    		continue;
    		
    	}
    	// Extract Groundstation and pull lat/long/height
    	
    	double[] lla = StationLookup.get(TDMData.get(0).get(0).station);
    	
    	GroundStation gst = null;
    	
    	// Currently no difference, although technically the coords of the stations are defined differently.
    	
    	if(TDMData.get(0).get(0).station.equals("SMART-01-A-SUTH") || TDMData.get(0).get(0).station.equals("SMART-01-B-SUTH"))
    	{
        	gst = new GroundStation(new TopocentricFrame(new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
    			    Constants.WGS84_EARTH_FLATTENING, FramesFactory.getITRF(ITRFVersion.ITRF_2000, IERSConventions.IERS_2010, false)),
        			new GeodeticPoint(lla[0]*Math.PI/180, lla[1]*Math.PI/180, lla[2]), TDMData.get(0).get(0).station));
    	}
    	else
    	{
        	gst = new GroundStation(new TopocentricFrame(new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
    			    Constants.WGS84_EARTH_FLATTENING, DataManager.getFrame("ITRF")),
        			new GeodeticPoint(lla[0]*Math.PI/180, lla[1]*Math.PI/180, lla[2]), TDMData.get(0).get(0).station));

    	}
    	
		gst.getPrimeMeridianOffsetDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
		gst.getPolarOffsetXDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
		gst.getPolarOffsetYDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
		
		if (referenceDate == null)
			referenceDate = DataManager.parseDateTime(TDMData.get(0).get(0).time);

    	// For each measurement time:
   	    	
    	    	
    	for(int i = 0; i < TDMData.size(); i++)
    	{
        	for(int j = 0; j < TDMData.get(i).size(); j++)
        	{
    	    	// Extract measurement/time

        		double RA = TDMData.get(i).get(j).rightAscension;
        		double Dec = TDMData.get(i).get(j).declination;
        		AbsoluteDate tMeas = DataManager.parseDateTime(TDMData.get(i).get(j).time);
        		
        		// convert UTC to gps week				
        		GNSSDate tMeasGPS = new GNSSDate(tMeas, SatelliteSystem.valueOf("GPS"));

        		int week = tMeasGPS.getWeekNumber();
        		int dayOfWeek = (int) tMeasGPS.getMilliInWeek()/(1000*60*60*24);
        		
        		// Read in GPS file corresponding to GPS week
        		
            	String SP3File = "D:\\School\\Work\\NATOData\\GNSS\\"+week+"\\igs"+week+dayOfWeek+".sp3";
            	
            	SP3File SP3ParsedData = null;
            	
            	try {
        			SP3ParsedData = new SP3Parser().parse(SP3File);
        			}
            	catch (IOException e) {
            		e.printStackTrace();
            		}
        		
        		// Pull Ephemeris for correct PRN
        	
            	SP3File.SP3Ephemeris SatelliteEphem = SP3ParsedData.getSatellites().get(PRN);

            	if(tMeas.compareTo(SatelliteEphem.getStop()) == 1)
            	{
            		System.out.println(week);
            		System.out.println(dayOfWeek);
            		System.out.println("Skipped Measurement!!");
            		continue;
            	}
            	
            	
            	
        		// Propagate to measurement time step
            	BoundedPropagator SP3Orbit = SatelliteEphem.getPropagator();
            	

            	SpacecraftState[] sta = new SpacecraftState[]{SP3Orbit.propagate(tMeas)};
            	

            	
        		// Produce measurement
            	
            	double[] obs = null;
            	double[] obsAzEl = {0,0};
            	/*
            	obs = new AngularRaDec(gst, FramesFactory.getEME2000(), tMeas, zeros, zeros, ones,
						   obssat).estimate(0, 0, sta).getEstimatedValue();
            	*/
            	
            	
            	if(TDMData.get(0).get(0).station.equals("ART"))
            	{
                	obs = new AngularRaDec(gst, FramesFactory.getGCRF(), tMeas, zeros, zeros, ones,
 						   obssat).estimate(0, 0, sta).getEstimatedValue();
            	}
            	else if(TDMData.get(0).get(0).station.equals("TRCK1") || TDMData.get(0).get(0).station.equals("TRCK2"))
            	{
                	obs = new AngularRaDec(gst, FramesFactory.getEME2000(), tMeas, zeros, zeros, ones,
 						   obssat).estimate(0, 0, sta).getEstimatedValue();
            	}
            	else if(TDMData.get(0).get(0).station.equals("SMART-01-A-SUTH") || TDMData.get(0).get(0).station.equals("SMART-01-B-SUTH"))
            	{
            		
                	double LTcorr = sta[0].getPVCoordinates().getPosition().subtract(gst.getBaseFrame().getPVCoordinates(tMeas, sta[0].getFrame()).getPosition()).getNorm() / Constants.SPEED_OF_LIGHT;
                	
                	sta = new SpacecraftState[]{SP3Orbit.propagate(new AbsoluteDate(tMeas, -LTcorr))};
    				
            		Vector3D relPosMinusLT = sta[0].getPVCoordinates(FramesFactory.getEME2000()).getPosition().subtract(gst.getBaseFrame().getPVCoordinates(tMeas, FramesFactory.getEME2000()).getPosition());
            				
    				Vector3D aberrationVel = gst.getBaseFrame().getPVCoordinates(tMeas, FramesFactory.getICRF()).getVelocity().add(gst.getBaseFrame().getPVCoordinates(tMeas, FramesFactory.getEME2000()).getVelocity());
            		
            		obs = new double[2];
            		obs[0] = relPosMinusLT.subtract(aberrationVel.scalarMultiply(LTcorr)).getAlpha();
            		obs[1] = relPosMinusLT.subtract(aberrationVel.scalarMultiply(LTcorr)).getDelta();
                	
            		
            	}
            	else
            	{
            		/*
                	obs = new AngularRaDec(gst, FramesFactory.getGCRF(), tMeas, zeros, zeros, ones,
 						   obssat).estimate(0, 0, sta).getEstimatedValue();
				   */
            		
                	
                	//Account for LT (very minor change in results...)
                	double LTcorr = sta[0].getPVCoordinates().getPosition().subtract(gst.getBaseFrame().getPVCoordinates(tMeas, sta[0].getFrame()).getPosition()).getNorm() / Constants.SPEED_OF_LIGHT;
                	
                	sta = new SpacecraftState[]{SP3Orbit.propagate(new AbsoluteDate(tMeas, -LTcorr))};
    				
            		Vector3D relPosMinusLT = sta[0].getPVCoordinates(FramesFactory.getGCRF()).getPosition().subtract(gst.getBaseFrame().getPVCoordinates(tMeas, FramesFactory.getGCRF()).getPosition());
            				
    				Vector3D aberrationVel = gst.getBaseFrame().getPVCoordinates(tMeas, FramesFactory.getICRF()).getVelocity().add(gst.getBaseFrame().getPVCoordinates(tMeas, FramesFactory.getGCRF()).getVelocity());
            		
            		obs = new double[2];
            		obs[0] = relPosMinusLT.subtract(aberrationVel.scalarMultiply(LTcorr)).getAlpha();
            		obs[1] = relPosMinusLT.subtract(aberrationVel.scalarMultiply(LTcorr)).getDelta();
                	
                	/*
            		AngularAzEl AzElModel = new AngularAzEl(gst, sta[0].getDate(), zeros, zeros, ones, obssat);
    				obsAzEl = AzElModel.estimate(0, 0, sta).getEstimatedValue();
            		GeodeticPoint gstPoint = gst.getOffsetGeodeticPoint(sta[0].getDate());

            		EarthStandardAtmosphereRefraction RefractionCorr = new EarthStandardAtmosphereRefraction();

            		double[] refractionRADEC = Conversion.convertAzElToRaDec(sta[0].getDate().toString(), 0, RefractionCorr.getRefraction(obsAzEl[1]), gstPoint.getLatitude(), gstPoint.getLongitude(), gstPoint.getAltitude(), "GCRF");
            		
            		obs[0] = obs[0] - refractionRADEC[0];
            		obs[1] = obs[1] - refractionRADEC[1];
            		*/
            		
            	}
            	

            	if(RA - obs[0] > Math.PI)
            	{
            		obs[0] = obs[0]+2*Math.PI;
            	}
            	else if(RA - obs[0] < -Math.PI)
            	{
            		obs[0] = obs[0]-2*Math.PI;
            	}

            	if(Dec - obs[1] > Math.PI)
            	{
            		obs[1] = obs[1]+2*Math.PI;
            	}
            	else if(Dec - obs[1] < -Math.PI)
            	{
            		obs[1] = obs[1]-2*Math.PI;
            	}
            	
            	// if file doesnt exist, output initial state
            	
            	if(!new File("D:\\School\\Work\\NATOData\\Output\\"+TDMData.get(0).get(0).NORAD+".txt").exists())
            	{
        			try {
        				
                    	String output = "D:\\School\\Work\\NATOData\\Output\\"+TDMData.get(0).get(0).NORAD+".txt";

                    	String data = sta[0].getPVCoordinates(FramesFactory.getEME2000()).toString()+"\n";
                    	
                    	Files.write(Paths.get(output), data.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.APPEND);
        				
        		
        			} catch (IOException e) {
        				
        					e.printStackTrace();
        			}
            	}
            	
            	// Save measurement/residual to file.
            	
        			try {
        				
                    	String output = "D:\\School\\Work\\NATOData\\Output\\"+TDMData.get(0).get(0).NORAD+".txt";

                    	String data = tMeas.toString()+","+TDMData.get(0).get(0).station+","+RA+","+Dec+","+obs[0]+","+obs[1]+"\n";
                    	
                    	Files.write(Paths.get(output), data.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.APPEND);
        				
        		
        			} catch (IOException e) {
        				
        					e.printStackTrace();
        			}

            	// Write JSON measurement to file.

    			try {
    				
                	String output = "D:\\School\\Work\\NATOData\\Output\\"+TDMData.get(0).get(0).NORAD+".json";

                	String data = "{\"Time\": \"" + tMeas.toString()+"Z\", \"Station\": \""+TDMData.get(0).get(0).station+"\", \"RightAscension\": "+RA+", \"Declination\": "+Dec + "},\n";
                	
                	Files.write(Paths.get(output), data.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.APPEND);
    				
    		
    			} catch (IOException e) {
    				
    					e.printStackTrace();
    			}
    			
    			// Write Files for each specific sensor.
        			
    			try {
    				
                	String output = "D:\\School\\Work\\NATOData\\Output\\StationData\\"+TDMData.get(0).get(0).station+".txt";

                	String data = tMeas.durationFrom(referenceDate)+","+RA+","+Dec+","+obs[0]+","+obs[1]+"\n";
                	
                	Files.write(Paths.get(output), data.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.APPEND);
    				
    		
    			} catch (IOException e) {
    				
    					e.printStackTrace();
    			}
    			
        	}
    		
    	}
    	
		}
    	
    }

    
   
    
}
