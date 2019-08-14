import os
import sys
import json
import math
import time
from astropy.time import Time
import numpy as np
import astrotools

def tocParse(rawData, fname, outdir, root):

    # Create sim_cfg.json file from raw obs_data.json file
    fnsim = fname + 'sim_cfg.json'
    # open up a reference sim_cfg file to edit
    with open(os.path.join(root, 'tle_sim_cfg.json'), "r") as fp:
        refSimFileData = json.load(fp)
    refSimFileData['Propagation']['Start'] = rawData['data']['obsList'][0]['time']
    refSimFileData['Propagation']['End'] = rawData['data']['obsList'][-1]['time']
    refSimFileData['Propagation']['InitialTLE'] = rawData['initialTle']
    refSimFileData['Stations'] = {
        rawData['sensor']['name'] : {'Latitude' : rawData['sensor']['location']['latitude']/180*math.pi,
                                    'Longitude' : rawData['sensor']['location']['longitude']/180*math.pi,
                                    'Altitude' : rawData['sensor']['location']['altitude']}
                                    }
    ## IF EXISTING FILTER OUTPUT FOR SENSOR SATELLITE COMBINATION
    idx = [i for i, char in enumerate(fname) if char == '_']
    sensorID = fname[idx[0]+1:idx[1]]
    satelliteID = fname[idx[1]+1:idx[2]]
    fnDir = os.listdir(outdir)
    if not fnDir:
        fnDir = ''
    else:
        fnValid = [fn for fn in fnDir if sensorID in fn and satelliteID in fn and 'od_out' in fn]
        if len(fnValid) != 0:
            with open(os.path.join(outdir, fnValid[-1]), "r") as fp:
                FilteredData = json.load(fp)
            refSimFileData['Propagation']['Start'] = FilteredData['Propagation']['Time']
    # Write sim_cfg file to dataTOC folder              
    with open(os.path.join(outdir,fnsim), 'w') as fp:
        json.dump(refSimFileData, fp, indent=4)

    # Create od_cfg.json file from raw obs_data.json file
    fnod = fname + 'od_cfg.json'
    # open up a reference od_cfg file to edit
    with open(os.path.join(root, 'tle_od_cfg.json'), "r") as fp:
        refOdFileData = json.load(fp)
    refOdFileData['Propagation']['Start'] = rawData['data']['obsList'][0]['time']
    refOdFileData['Propagation']['End'] = rawData['data']['obsList'][-1]['time']
    refOdFileData['Propagation']['InitialTLE'] = rawData['initialTle']
    refOdFileData['Stations'] = {
        rawData['sensor']['name'] : {'Latitude' : rawData['sensor']['location']['latitude']/180*math.pi,
                                    'Longitude' : rawData['sensor']['location']['longitude']/180*math.pi,
                                    'Altitude' : rawData['sensor']['location']['altitude']}
                                    }
    # ## IF EXISTING FILTER OUTPUT FOR SENSOR SATELLITE COMBINATION
    if not fnDir:
        fnDir = ''
    else:
        fnValid = [fn for fn in fnDir if sensorID in fn and satelliteID in fn and 'od_out' in fn]
        if len(fnValid) != 0:
            refOdFileData['Propagation'] = {'Start' : FilteredData['Propagation']['Time'],  
                                            'End' : rawData['data']['obsList'][-1]['time'],
                                            'InitialState' : FilteredData['Propagation']['State'][:6]}
            cov = np.diag(np.array(FilteredData['Estimation'][-1]['EstimatedCovariance']))
            refOdFileData['Estimation']['Covariance'] = cov.tolist()
    # Write od_cfg file to dataTOC folder 
    with open(os.path.join(outdir,fnod), 'w') as fp:
        json.dump(refOdFileData, fp, indent=4)

    return

# Get orbdetpy folder path
odpdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Get path to data folder
datdir = os.path.join(odpdir, "examples", "dataTOC")
# Append orbdetpy path to call java methods
sys.path.append(odpdir)
from orbdetpy import determineOrbit, simulateMeasurements

print("run_TOC start batch: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
# Parse raw obs_data.json in dataTOC folder
for root, dirs, files in os.walk(datdir):
    outdir = os.path.join(root, "output")
    if 'output' in root:
        break
    # Create a sim_cfg.json and od_cfg.json file for raw obs_data.json file
    for fname in files:
        idx = fname.find("obs_data.json")
        if (idx == -1):
            continue
        fname = fname[:idx]
        print("Processing {}".format(fname))
        print("%s" % time.strftime("%Y-%m-%d %H:%M:%S"))
        with open(os.path.join(root, fname+"obs_data.json"), "r") as fp:
            rawData = json.load(fp)
        tocParse(rawData, fname, outdir, root)
        # Create obs.json files of simulated measurement and reference state
        with open(os.path.join(outdir, fname+"sim_cfg.json"), "r") as fp:
            config = fp.read()
        obs = simulateMeasurements(config)
        with open(os.path.join(outdir, fname + "obs.json"), "w") as fp:
            fp.write(obs)
        # Interpolate simulated state and create a obs_data.json file with the real obs and times
        with open(os.path.join(outdir, fname + "obs.json"), "r") as fp:
            simobs = json.load(fp)
        with open(os.path.join(root, fname + "obs_data.json"), "r") as fp:
            realobs = json.load(fp)
        output = []
        for real in realobs["data"]["obsList"]:
            output.append({"TrueState" : {"Cartesian" : [0,0,0,0,0,0]}, "Time" : real['time'], "Station" : real['station'],
                            "RightAscension" : real['rightAscension'], "Declination" : real['declination']})
        pvRef = [est['TrueState']['Cartesian'] for est in simobs]
        timesim = [est['Time'] for est in simobs]
        timesim = Time(timesim, format='isot', scale='utc')
        timesim = timesim.mjd
        newtime = [est['time'] for est in realobs['data']['obsList']]
        newtime = Time(newtime, format='isot', scale='utc')
        newtime = newtime.mjd
        degree = 10
        pvObs = astrotools.interpolationPV(newtime, timesim, np.array(pvRef), degree)
        count = 0
        for state in output:
            state['TrueState']['Cartesian'] = pvObs[count,:].tolist()
            count += 1
        with open(os.path.join(outdir, fname + "obs_data.json"), 'w') as fp:
            json.dump(output, fp, indent=4)
        # Run orbit determination and create od_output.json files
        with open(os.path.join(outdir, fname + "obs_data.json"), "r") as fp:
            obs = fp.read()
        with open(os.path.join(outdir, fname + "od_cfg.json"), "r") as fp:
            config = json.load(fp)
        fit = determineOrbit(json.dumps(config), obs)
        with open(os.path.join(outdir, fname + "od_out.json"), "w") as fp:
            fp.write(fit)
        
        output = []
        exact = [x for x in json.loads(obs) if ("Station" in x or "PositionVelocity" in x)]
        estim = json.loads(fit)["Estimation"]
        for exa, est in zip(exact, estim):
            if ("TrueState" in exa):
                X = exa["TrueState"]["Cartesian"][:6]
            else:
                X = exa["PositionVelocity"][:6]
            diff = [X[i] - est["EstimatedState"][i] for i in range(6)]
            sigm = [3.0*math.sqrt(est["EstimatedCovariance"][i][i]) for i in range(6)]
            output.append({"Time" : exa["Time"], "Station" : exa.get("Station"),
                            "StateResidual" : diff, "Covariance3Sigma" : sigm,
                            "WithinBounds" : [diff[i] > -sigm[i] and diff[i] < sigm[i] for i in range(6)]})
        with open(os.path.join(outdir, fname+"_diff.json"), "w") as fp:
            json.dump(output, fp, indent = 4)
        print("%s" % time.strftime("%Y-%m-%d %H:%M:%S"))
        os.remove(os.path.join(root, fname + "obs_data.json"))
        os.remove(os.path.join(outdir, fname + "sim_cfg.json")) 
        os.remove(os.path.join(outdir, fname + "obs.json"))
        
print("run_TOC end batch  : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))