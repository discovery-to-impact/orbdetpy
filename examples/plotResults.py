import plotly
import plotly.graph_objects as go
import json
import numpy as np
from plotly.subplots import make_subplots
import dateutil.parser
import math
import scipy.io as sio
from scipy.interpolate import interp1d
import astrotools
import plotly.figure_factory as ff
import time as tsys
import plotly.io as pio
pio.templates.default = "none"

def plot(obs_data, od_cfg, od_out, smootherRun,plots):
    """Plots results from OD
    :param mT: str before _obs,_od json files
    :param smootherRun: 1 or 0 to run smoother 
    :param plots: 1 or 0 to plot (only create table)
    """
    print("Plot start : %s" % tsys.strftime("%Y-%m-%d %H:%M:%S"))

    plotRef = 1

    mT = obs_data[5:9]

    if smootherRun:
        print("Smoother start : %s" % tsys.strftime("%Y-%m-%d %H:%M:%S"))
        xPre, xPost, PPre, PPost, SPpre = astrotools.readODoutJSON(od_out)
        xS, PS, FSCT, time = astrotools.rtsSmoother(xPre, xPost, PPre, PPost, SPpre)
        FSCT = abs(FSCT)
        sigPs = astrotools.get3sig(PS)
        print("Smoother end : %s" % tsys.strftime("%Y-%m-%d %H:%M:%S"))
   
    with open(od_cfg, "r") as fp:
        cfg = json.load(fp)
    with open(obs_data, "r") as fp:
        inp = json.load(fp)
    with open(od_out, "r") as fp:
        out = json.load(fp)["Estimation"]

    key = tuple(cfg["Measurements"].keys())
    dmcrun = (cfg["Estimation"].get("DMCCorrTime", 0.0) > 0.0 and
                cfg["Estimation"].get("DMCSigmaPert", 0.0) > 0.0)

    tstamp, prefit, posfit, inocov, params, estmacc, estmcov = [], [], [], [], [], [], []
    pvRef, pvEst, tstr, ang1, ang2, estcov, sigestcov, pvPre, precov, preSP = [], [], [], [], [], [], [], [], [], []
    
    # Get site location
    key2 = tuple(cfg["Stations"].keys())
    lat = cfg['Stations'][key2[0]]['Latitude']
    lon = cfg['Stations'][key2[0]]['Longitude']
    alt = cfg['Stations'][key2[0]]['Altitude']

    for i, o in zip(inp, out):
        tstamp.append(dateutil.parser.isoparse(i["Time"]))

        # Get time string
        tstr.append(i["Time"])
        # Get reference state
        PV = np.array(i["TrueState"]["Cartesian"][:6])
        pvRef.append(PV)
        # Get filter state
        PV = np.array(o["EstimatedState"][:6])
        pvEst.append(PV)
        # Get filter covariance
        covi = np.array(o["EstimatedCovariance"])
        estcov.append(covi)
        # Get filter covariance 3-sigma diagonal
        p = []
        for m in range(len(o["EstimatedCovariance"])):
            p.append(3.0*np.sqrt(o["EstimatedCovariance"][m][m]))
        sigestcov.append(p)
        # Get smoother inputs
        if smootherRun:
            # Get pre filter state
            PV = np.array(o["xBar"][:6])
            pvPre.append(PV) 
            # Get pre filter covariance
            covi = np.array(o["PBar"])
            precov.append(covi)
            # Get pre filter sigma points
            covi = np.array(o["SPpre"])
            preSP.append(covi)

        if ("PositionVelocity" in key):
            prefit.append([ix - ox for ix, ox in zip(i["PositionVelocity"],
                                                        o["PreFit"]["PositionVelocity"])])
            posfit.append([ix - ox for ix, ox in zip(i["PositionVelocity"],
                                                        o["PostFit"]["PositionVelocity"])])
        else:
            # Get obs
            ang1.append(i[key[0]])
            ang2.append(i[key[1]])
            prefit.append([i[key[0]] - o["PreFit"][key[0]][-1],
                            i[key[1]] - o["PreFit"][key[1]][-1]])
            posfit.append([i[key[0]] - o["PostFit"][key[0]][-1],
                            i[key[1]] - o["PostFit"][key[1]][-1]])

        p = []
        for m in range(len(o["InnovationCovariance"])):
            p.append(3.0*np.sqrt(o["InnovationCovariance"][m][m]))
        inocov.append(p)

        if (len(o["EstimatedState"]) > 6):
            if (dmcrun):
                params.append(o["EstimatedState"][6:-3])
            else:
                params.append(o["EstimatedState"][6:])

        if (dmcrun):
            r = np.array(o["EstimatedState"][:3])
            r /= np.linalg.norm(r)
            v = np.array(o["EstimatedState"][3:6])
            v /= np.linalg.norm(v)
            h = np.cross(r, v)
            rot = np.vstack((r, np.cross(h, r), h))
            estmacc.append(rot.dot(o["EstimatedState"][-3:]))

    pre = np.array(prefit)
    pos = np.array(posfit)
    cov = np.array(inocov)
    par = np.array(params)
    estmacc = np.array(estmacc)
    tim = [(t - tstamp[0]).total_seconds()/3600 for t in tstamp]

    obs = np.array([ang1,ang2])
    refobs = np.array(pvRef)

    xPost = np.array(pvEst)
    PPost = np.array(estcov)
    sigPf = np.array(sigestcov)
    sigPf = sigPf.T

    xPre = np.array(pvPre)
    PPre = np.array(precov)
    SPpre = np.array(preSP)

    numObs = np.size(xPost[:,0])

    angles = ("Azimuth", "Elevation", "RightAscension", "Declination")
    if (key[0] in angles and key[1] in angles):
        pre *= 648000.0/math.pi
        pos *= 648000.0/math.pi
        cov *= 648000.0/math.pi
        units = ("[arcsec]", "[arcsec]")
    else:
        if ("PositionVelocity" in key):
            units = ("[m]", "[m]", "[m]", "[m/s]", "[m/s]", "[m/s]")
        else:
            units = ("[m]", "[m/s]")

    if ("PositionVelocity" in key):
        ylabs = (r"$\Delta x$", r"$\Delta y$", r"$\Delta z$",
                    r"$\Delta v_x$", r"$\Delta v_y$", r"$\Delta v_z$")
        order = (1, 3, 5, 2, 4, 6)
    else:
        ylabs = key
    #############################################################################################

    # Get estimated measurement of smoother states and reference states
    smoothMeas, smoothRes, refMeas, refRes = [], [], [], []
    sig1 = np.array(cfg['Measurements'][key[0]]['Error'])
    sig2 = np.array(cfg['Measurements'][key[1]]['Error'])
    for ii in range(0,numObs,1):
        ti = tstr[ii]
        sigma = [sig1,sig2]
        if plotRef:
            pv = refobs[ii,:].tolist()
            angular = obs[:,ii].tolist()
            if key[0] == 'RightAscension':
                output = astrotools.pv2radec(od_cfg, lat, lon, alt, ti, angular, sigma, pv)
            if key[0] == 'Azimuth':
                output = astrotools.pv2azel(od_cfg, lat, lon, alt, ti, angular, sigma, pv)
            if key[0] == 'Range':
                range1 = obs[0,ii].tolist()
                rangerate1 = obs[1,ii].tolist()
                output1 = astrotools.pv2range(od_cfg, lat, lon, alt, ti, range1, sig1, pv)
                output2 = astrotools.pv2rangerate(od_cfg, lat, lon, alt, ti, rangerate1, sig2, pv)
                output = [output1,output2]
            refMeas.append(output)
            output = np.squeeze(np.array(output))
            angular = np.array(angular)
            refRes.append(angular - output)
        
        if smootherRun:
            pv = xS[0:6,ii].tolist()
            angular = obs[:,ii].tolist()
            if key[0] == 'RightAscension':
                output = astrotools.pv2radec(od_cfg, lat, lon, alt, ti, angular, sigma, pv)
            if key[0] == 'Azimuth':
                output = astrotools.pv2azel(od_cfg, lat, lon, alt, ti, angular, sigma, pv)
            if key[0] == 'Range':
                range1 = obs[0,ii].tolist()
                rangerate1 = obs[1,ii].tolist()
                output1 = astrotools.pv2range(od_cfg, lat, lon, alt, ti, range1, sig1, pv)
                output2 = astrotools.pv2rangerate(od_cfg, lat, lon, alt, ti, rangerate1, sig2, pv)
                output = [output1,output2]
            smoothMeas.append(output)
            output = np.squeeze(np.array(output))
            angular = np.array(angular)
            smoothRes.append(angular - output)

    refMeas = np.array(refMeas)
    refRes = np.array(refRes)
    smoothMeas = np.array(smoothMeas)
    smoothRes = np.array(smoothRes)
    if key[0] in ('RightAscension','Azimuth'):
        refRes = refRes*648000.0/math.pi
        smoothRes = smoothRes*648000.0/math.pi
        sig1 = sig1*648000.0/math.pi
        sig2 = sig2*648000.0/math.pi
    # Get RIC from cartesian reference positions
    xrefRIC, vrefRIC = [], []
    for ii in range(0,numObs,1):
        xrefRIC.append(astrotools.xECI2RIC(refobs[ii,0:3].T, refobs[ii,0:3].T, refobs[ii,3:6].T))
        vrefRIC.append(astrotools.xECI2RIC(refobs[ii,3:6].T, refobs[ii,0:3].T, refobs[ii,3:6].T))
    xrefRIC = np.array(xrefRIC)
    vrefRIC = np.array(vrefRIC)
    # Get RIC from cartesian filter positions and covariance
    xRIC, vRIC, PRIC = [], [], []
    for ii in range(0,numObs,1):
        xRIC.append(astrotools.xECI2RIC(xPost[ii,0:3].T, refobs[ii,0:3].T, refobs[ii,3:6].T))
        vRIC.append(astrotools.xECI2RIC(xPost[ii,3:6].T, refobs[ii,0:3].T, refobs[ii,3:6].T))
        PRIC.append(astrotools.PECI2RIC(PPost[ii,0:6,0:6], refobs[ii,0:3].T, refobs[ii,3:6].T))
    xRIC = np.array(xRIC)
    vRIC = np.array(vRIC)
    PRIC = np.array(PRIC)
    PRIC = PRIC.T
    sigPRIC = astrotools.get3sig(PRIC)
    # Get RIC from cartesian smoother positions and covariance
    if smootherRun:
        xSRIC, vSRIC, PSRIC = [], [], []
        for ii in range(0,numObs,1):
            xSRIC.append(astrotools.xECI2RIC(xS[0:3,ii], refobs[ii,0:3].T, refobs[ii,3:6].T))
            vSRIC.append(astrotools.xECI2RIC(xS[3:6,ii], refobs[ii,0:3].T, refobs[ii,3:6].T))
            PSRIC.append(astrotools.PECI2RIC(PS[0:6,0:6,ii], refobs[ii,0:3].T, refobs[ii,3:6].T))
        xSRIC = np.array(xSRIC)
        vSRIC = np.array(vSRIC)
        PSRIC = np.array(PSRIC)
        PSRIC = PSRIC.T
        sigPSRIC = astrotools.get3sig(PSRIC)

    # STATE RESIDUALS (REFERNCE STATE - FILTER STATE)
    PVDiff = refobs[:,0:6] - xPost[:,0:6]
    meanstdPVfilter = astrotools.meanstdPV(PVDiff)
    # STATE RESIDUALS (REFERENCE STATE - SMOOTHER STATE)
    if smootherRun:
        posvelDiff = refobs[:,0:6] - xS[0:6,:].T
        meanstdPVsmooth = astrotools.meanstdPV(posvelDiff)
    # STATE RESIDUALS (RIC REFERNCE STATE - RIC FILTER STATE)
    xDiffRIC = np.zeros([numObs,6])
    xDiffRIC[:,0:3] = xrefRIC - xRIC
    xDiffRIC[:,3:6] = vrefRIC - vRIC
    meanstdRICfilter = astrotools.meanstdPV(xDiffRIC)
    # STATE RESIDUALS (RIC REFERNCE STATE - RIC SMOOTHER STATE)
    if smootherRun:
        xSDiffRIC = np.zeros([numObs,6])
        xSDiffRIC[:,0:3] = xrefRIC - xSRIC
        xSDiffRIC[:,3:6] = vrefRIC - vSRIC
        meanstdRICsmooth = astrotools.meanstdPV(xSDiffRIC)
    # MEAN & STD (REFERENCE MEASUREMENT RESIDUALS)
    if plotRef:
        meanstdANGref = astrotools.meanstdANG(refRes)
    # MEAN & STD (FILTER MEASUREMENT RESIDUALS)
    meanstdANGpre = astrotools.meanstdANG(pre)
    meanstdANGpost = astrotools.meanstdANG(pos)
    # MEAN & STD (SMOOTHER MEASUREMENT RESIDUALS)
    if smootherRun:
        meanstdANGsmooth = astrotools.meanstdANG(smoothRes)

    print("Plot end : %s" % tsys.strftime("%Y-%m-%d %H:%M:%S"))
    
    ################################################## PLOTLY #############################################################

    initialCov = np.array(cfg['Estimation']['Covariance'])
    initialProcNoise = np.array(cfg['Estimation']['ProcessNoise'])
    initialANGNoise = [sig1,sig2]
    
    data_matrix = [['DESCRIPTION',' ']]
    data_matrix.append(['Name', mT])
    data_matrix.append(['Filter Type', cfg['Estimation']['Filter']])
    data_matrix.append(['Number of observations', numObs])
    data_matrix.append(['Propagation Start time [UTC]', cfg['Propagation']['Start']])
    data_matrix.append(['Propagation End time [UTC]', cfg['Propagation']['End']])
    data_matrix.append(['Initial Covariance [m,m/s]', initialCov])
    data_matrix.append(['Initial Process Noise [m,m/s]', initialProcNoise])
    data_matrix.append(['Initial Measurement Noise '+key[0]+' '+units[0], initialANGNoise[0]])
    data_matrix.append(['Initial Measurement Noise '+key[1]+' '+units[1], initialANGNoise[1]])
    
    data_matrix.append(['PLOTS', ' '])
    data_matrix.append(['Filter Measurement Residuals '+key[0]+' '+units[0], 'Prefit: Mean = ' + format(meanstdANGpre[0],'.2f') + ' | STD = ' + format(meanstdANGpre[1],'.2f')])
    data_matrix.append(['Filter Measurement Residuals '+key[1]+' '+units[1], 'Prefit: Mean = ' + format(meanstdANGpre[2],'.2f') + ' | STD = ' + format(meanstdANGpre[3],'.2f')])
    data_matrix.append(['Filter Measurement Residuals '+key[0]+' '+units[0], 'Postfit: Mean = ' + format(meanstdANGpost[0],'.2f') + ' | STD = ' + format(meanstdANGpost[1],'.2f')])
    data_matrix.append(['Filter Measurement Residuals '+key[1]+' '+units[1], 'Postfit: Mean = ' + format(meanstdANGpost[2],'.2f') + ' | STD = ' + format(meanstdANGpost[3],'.2f')])
    if plotRef:
        data_matrix.append(['Filter Measurement Residuals '+key[0]+' '+units[0], 'Reference: Mean = ' + format(meanstdANGref[0],'.2f') + ' | STD = ' + format(meanstdANGref[1],'.2f')])
        data_matrix.append(['Filter Measurement Residuals '+key[1]+' '+units[1], 'Reference: Mean = ' + format(meanstdANGref[2],'.2f') + ' | STD = ' + format(meanstdANGref[3],'.2f')])
   
    if smootherRun:
        data_matrix.append(['Smoother Measurements Residuals '+key[0]+' '+units[0], 'Smoothed: Mean = ' + format(meanstdANGsmooth[0],'.2f') + ' | STD = ' + format(meanstdANGsmooth[1],'.2f')])
        data_matrix.append(['Smoother Measurements Residuals '+key[1]+' '+units[1], 'Smoothed: Mean = ' + format(meanstdANGsmooth[2],'.2f') + ' | STD = ' + format(meanstdANGsmooth[3],'.2f')])
    
    data_matrix.append(['Reference - Filter Position Residuals [m]', 'Position x: Mean = ' + format(meanstdPVfilter[0],'.2f') + ' | STD = ' + format(meanstdPVfilter[1],'.2f')])
    data_matrix.append(['Reference - Filter Position Residuals [m]', 'Position y: Mean = ' + format(meanstdPVfilter[2],'.2f') + ' | STD = ' + format(meanstdPVfilter[3],'.2f')])
    data_matrix.append(['Reference - Filter Position Residuals [m]', 'Position z: Mean = ' + format(meanstdPVfilter[4],'.2f') + ' | STD = ' + format(meanstdPVfilter[5],'.2f')])
    # data_matrix.append(['Reference - Filter Position and Velocity Residuals [m/s]', 'Position vx: Mean = ' + format(meanstdPVfilter[6],'.2f') + ' | STD = ' + format(meanstdPVfilter[7],'.2f')])
    # data_matrix.append(['Reference - Filter Position and Velocity Residuals [m/s]', 'Position vy: Mean = ' + format(meanstdPVfilter[8],'.2f') + ' | STD = ' + format(meanstdPVfilter[9],'.2f')])
    # data_matrix.append(['Reference - Filter Position and Velocity Residuals [m/s]', 'Position vz: Mean = ' + format(meanstdPVfilter[10],'.2f') + ' | STD = ' + format(meanstdPVfilter[11],'.2f')])
    if smootherRun:
        data_matrix.append(['Reference - Smoother Position Residuals [m]', 'Position x: Mean = ' + format(meanstdPVsmooth[0],'.2f') + ' | STD = ' + format(meanstdPVsmooth[1],'.2f')])
        data_matrix.append(['Reference - Smoother Position Residuals [m]', 'Position y: Mean = ' + format(meanstdPVsmooth[2],'.2f') + ' | STD = ' + format(meanstdPVsmooth[3],'.2f')])
        data_matrix.append(['Reference - Smoother Position Residuals [m]', 'Position z: Mean = ' + format(meanstdPVsmooth[4],'.2f') + ' | STD = ' + format(meanstdPVsmooth[5],'.2f')])
    # data_matrix.append(['Reference - Smoother Position and Velocity Residuals [m/s]', 'Position vx: Mean = ' + format(meanstdPVsmooth[6],'.2f') + ' | STD = ' + format(meanstdPVsmooth[7],'.2f')])
    # data_matrix.append(['Reference - Smoother Position and Velocity Residuals [m/s]', 'Position vy: Mean = ' + format(meanstdPVsmooth[8],'.2f') + ' | STD = ' + format(meanstdPVsmooth[9],'.2f')])
    # data_matrix.append(['Reference - Smoother Position and Velocity Residuals [m/s]', 'Position vz: Mean = ' + format(meanstdPVsmooth[10],'.2f') + ' | STD = ' + format(meanstdPVsmooth[11],'.2f')])
    
    data_matrix.append(['RIC Reference - Filter Position Residuals [m]', 'Position Radial: Mean = ' + format(meanstdRICfilter[0],'.2f') + ' | STD = ' + format(meanstdRICfilter[1],'.2f')])
    data_matrix.append(['RIC Reference - Filter Position Residuals [m]', 'Position In-track: Mean = ' + format(meanstdRICfilter[2],'.2f') + ' | STD = ' + format(meanstdRICfilter[3],'.2f')])
    data_matrix.append(['RIC Reference - Filter Position Residuals [m]', 'Position Cross-track: Mean = ' + format(meanstdRICfilter[4],'.2f') + ' | STD = ' + format(meanstdRICfilter[5],'.2f')])
    if smootherRun:
        data_matrix.append(['RIC Reference - Smoother Position Residuals [m]', 'Position Radial: Mean = ' + format(meanstdRICsmooth[0],'.2f') + ' | STD = ' + format(meanstdRICsmooth[1],'.2f')])
        data_matrix.append(['RIC Reference - Smoother Position Residuals [m]', 'Position In-track: Mean = ' + format(meanstdRICsmooth[2],'.2f') + ' | STD = ' + format(meanstdRICsmooth[3],'.2f')])
        data_matrix.append(['RIC Reference - Smoother Position Residuals [m]', 'Position Cross-track: Mean = ' + format(meanstdRICsmooth[4],'.2f') + ' | STD = ' + format(meanstdRICsmooth[5],'.2f')])

    table = ff.create_table(data_matrix, index=True)

    # Make text size larger
    for i in range(len(table.layout.annotations)):
        table.layout.annotations[i].font.size = 20

    # Execute
    plotly.offline.plot(table, filename='Plots/index_table.html', auto_open=True)

    ########################################################################################
    if plots:
        TITLE = 'Filter Measurement Residuals'
        # Plot Filter Measurement Residuals
        title1 = 'Prefit: Mean = ' + format(meanstdANGpre[0],'.2f') + ' | STD = ' + format(meanstdANGpre[1],'.2f')
        trace1 = go.Scatter(x=tstr, y=pre[:,0], mode = 'markers', name = key[0], marker = dict(color = 'rgb(51, 153, 255)'))
        title2 = 'Prefit: Mean = ' + format(meanstdANGpre[2],'.2f') + ' | STD = ' + format(meanstdANGpre[3],'.2f')
        trace2 = go.Scatter(x=tstr, y=pre[:,1], mode = 'markers', name = key[1], marker = dict(color = 'rgb(0, 204, 0)'))
        title3 = 'Postfit: Mean = ' + format(meanstdANGpost[0],'.2f') + ' | STD = ' + format(meanstdANGpost[1],'.2f')
        trace3 = go.Scatter(x=tstr, y=pos[:,0], mode = 'markers', name = key[0], marker = dict(color = 'rgb(51, 153, 255)'), showlegend = False)
        title4 = 'Postfit: Mean = ' + format(meanstdANGpost[2],'.2f') + ' | STD = ' + format(meanstdANGpost[3],'.2f')
        trace4 = go.Scatter(x=tstr, y=pos[:,1], mode = 'markers', name = key[1], marker = dict(color = 'rgb(0, 204, 0)'), showlegend = False)
        # Plot 3-sigma Bounds
        trace5 = go.Scatter(x=tstr, y=cov[:,0], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds')
        trace6 = go.Scatter(x=tstr, y=-cov[:,0], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace7 = go.Scatter(x=tstr, y=cov[:,1], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace8 = go.Scatter(x=tstr, y=-cov[:,1], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        if plotRef:
            # Plot Reference residuals
            title5 = 'Reference: Mean = ' + format(meanstdANGref[0],'.2f') + ' | STD = ' + format(meanstdANGref[1],'.2f')
            trace9 = go.Scatter(x=tstr, y=refRes[:,0], mode = 'markers', name = key[0], marker = dict(color = 'rgb(51, 153, 255)'), showlegend = False)
            title6 = 'Reference: Mean = ' + format(meanstdANGref[2],'.2f') + ' | STD = ' + format(meanstdANGref[3],'.2f')
            trace10 = go.Scatter(x=tstr, y=refRes[:,1], mode = 'markers', name = key[1], marker = dict(color = 'rgb(0, 204, 0)'), showlegend = False)
        # Create subplot
        fig = make_subplots(rows=2, cols=2, subplot_titles=(title1, title2, title3, title4), shared_xaxes=True)
        if plotRef:
            fig = make_subplots(rows=3, cols=2, subplot_titles=(title1, title2, title3, title4, title5, title6), shared_xaxes=True)
        # Assign traces to subplots
        fig.append_trace(trace1, 1, 1)
        fig.append_trace(trace2, 1, 2)
        fig.append_trace(trace3, 2, 1)
        fig.append_trace(trace4, 2, 2)
        fig.append_trace(trace5, 2, 1)
        fig.append_trace(trace6, 2, 1)
        fig.append_trace(trace7, 2, 2)
        fig.append_trace(trace8, 2, 2)
        if plotRef:
            fig.append_trace(trace9, 3, 1)
            fig.append_trace(trace10, 3, 2)
        # Configure Layout
        fig['layout'].update(title= TITLE)
        fig['layout']['xaxis1'].update(title='Time since epoch [hrs]')
        fig['layout']['xaxis2'].update(title='Time since epoch [hrs]')
        fig['layout']['yaxis1'].update(title=key[0]+' '+units[0])
        fig['layout']['yaxis2'].update(title=key[1]+' '+units[1])
        fig['layout']['yaxis3'].update(title=key[0]+' '+units[0], range = [-5,5])
        fig['layout']['yaxis4'].update(title=key[1]+' '+units[1], range = [-5,5])
        if plotRef:
            fig['layout']['yaxis5'].update(title=key[0]+' '+units[0])
            fig['layout']['yaxis6'].update(title=key[1]+' '+units[1])
        # Execute
        plotly.offline.plot(fig, filename='Plots/filtermeas.html', auto_open=True)
        # plotly.offline.plot(fig, image = 'png', image_filename='filtermeas', output_type='file', image_width=1920, image_height=1080 , validate=False)

        ########################################################################################

        if smootherRun:
            TITLE = 'Smoother Measurements Residuals'
            # Plot Smoother residuals
            title1 = 'Smoothed: Mean = ' + format(meanstdANGsmooth[0],'.2f') + ' | STD = ' + format(meanstdANGsmooth[1],'.2f')
            trace1 = go.Scatter(x=tstr, y=smoothRes[:,0], mode = 'markers', name = 'residuals', marker = dict(color = 'rgb(51, 153, 255)'))
            title2 = 'Smoothed: Mean = ' + format(meanstdANGsmooth[2],'.2f') + ' | STD = ' + format(meanstdANGsmooth[3],'.2f')
            trace2 = go.Scatter(x=tstr, y=smoothRes[:,1], mode = 'markers', name = 'residuals', marker = dict(color = 'rgb(0, 204, 0)'))
            # Plot 3-sigma Bounds
            # trace3 = go.Scatter(x=tim, y=cov[:,0], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds')
            # trace4 = go.Scatter(x=tim, y=-cov[:,0], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            # trace5 = go.Scatter(x=tim, y=cov[:,1], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            # trace6 = go.Scatter(x=tim, y=-cov[:,1], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            # Plot FSCT
            title3 = 'McReynolds Smoother Filter Consistency Check: Position [m]'
            trace7 = go.Scatter(x=tstr, y=FSCT[0,:], mode = 'markers', name = 'x', marker = dict(color = 'rgb(51, 153, 255)'))
            trace8 = go.Scatter(x=tstr, y=FSCT[1,:], mode = 'markers', name = 'y', marker = dict(color = 'rgb(255, 192, 56)'))
            trace9 = go.Scatter(x=tstr, y=FSCT[2,:], mode = 'markers', name = 'z', marker = dict(color = 'rgb(0, 204, 0)'))
            title4 = 'McReynolds Smoother Filter Consistency Check: Velocity [m/s]'
            trace10 = go.Scatter(x=tstr, y=FSCT[3,:], mode = 'markers', name = 'vx', marker = dict(color = 'rgb(51, 153, 255)'))
            trace11 = go.Scatter(x=tstr, y=FSCT[4,:], mode = 'markers', name = 'vy', marker = dict(color = 'rgb(255, 192, 56)'))
            trace12 = go.Scatter(x=tstr, y=FSCT[5,:], mode = 'markers', name = 'vz', marker = dict(color = 'rgb(0, 204, 0)'))
            # Create subplot
            fig = make_subplots(rows=2, cols=2, subplot_titles=(title1, title2, title3, title4), shared_xaxes=True)
            # Assign traces to subplots
            fig.append_trace(trace1, 1, 1)
            fig.append_trace(trace2, 1, 2)
            # fig.append_trace(trace3, 1, 1)
            # fig.append_trace(trace4, 1, 1)
            # fig.append_trace(trace5, 1, 2)
            # fig.append_trace(trace6, 1, 2)
            fig.append_trace(trace7, 2, 1)
            fig.append_trace(trace8, 2, 1)
            fig.append_trace(trace9, 2, 1)
            fig.append_trace(trace10, 2, 2)
            fig.append_trace(trace11, 2, 2)
            fig.append_trace(trace12, 2, 2)
            # Configure Layout
            fig['layout'].update(title= TITLE)
            fig['layout']['xaxis1'].update(title='Time since epoch [hrs]')
            fig['layout']['xaxis2'].update(title='Time since epoch [hrs]')
            fig['layout']['yaxis1'].update(title=key[0]+' '+units[0])
            fig['layout']['yaxis2'].update(title=key[1]+' '+units[1])
            fig['layout']['yaxis3'].update(title='FSCT')
            fig['layout']['yaxis4'].update(title='FSCT')
            # Execute
            plotly.offline.plot(fig, filename='Plots/smoothermeas.html', auto_open=True)
            # plotly.offline.plot(fig, image = 'png', image_filename='smoothermeas', output_type='file', image_width=1920, image_height=1080 , validate=False)

        ########################################################################################
        if smootherRun:
            TITLE = 'Estimated Parameters'
            # Plot Smoother residuals
            title1 = 'Filter Coefficient of Drag'
            trace1 = go.Scatter(x=tstr, y=par[:,0], mode = 'markers', name = 'Cd', marker = dict(color = 'rgb(51, 153, 255)'))
            title2 = 'Filter Coefficient of Reflectivity'
            trace2 = go.Scatter(x=tstr, y=par[:,1], mode = 'markers', name = 'Cr', marker = dict(color = 'rgb(0, 204, 0)'))
            title3 = 'Smoother Coefficient of Drag'
            trace3 = go.Scatter(x=tstr, y=xS[6,:], mode = 'markers', name = 'Cd', marker = dict(color = 'rgb(51, 153, 255)'))
            title4 = 'Smoother Coefficient of Reflectivity'
            trace4 = go.Scatter(x=tstr, y=xS[7,:], mode = 'markers', name = 'Cr', marker = dict(color = 'rgb(0, 204, 0)'))
            # Plot FSCT
            title5 = 'McReynolds Smoother Filter Consistency Check: Coefficient of Drag'
            trace5 = go.Scatter(x=tstr, y=FSCT[6,:], mode = 'markers', name = 'Cd', marker = dict(color = 'rgb(51, 153, 255)'), showlegend = False)
            title6 = 'McReynolds Smoother Filter Consistency Check: Coefficient of Reflectivity'
            trace6 = go.Scatter(x=tstr, y=FSCT[7,:], mode = 'markers', name = 'Cr', marker = dict(color = 'rgb(0, 204, 0)'), showlegend = False)
            # Create subplot
            fig = make_subplots(rows=3, cols=2, subplot_titles=(title1, title2, title3, title4, title5, title6), shared_xaxes=True)
            # Assign traces to subplots
            fig.append_trace(trace1, 1, 1)
            fig.append_trace(trace2, 1, 2)
            fig.append_trace(trace3, 2, 1)
            fig.append_trace(trace4, 2, 2)
            fig.append_trace(trace5, 3, 1)
            fig.append_trace(trace6, 3, 2)
            # Configure Layout
            fig['layout'].update(title= TITLE)
            fig['layout']['xaxis1'].update(title='Time since epoch [hrs]')
            fig['layout']['xaxis2'].update(title='Time since epoch [hrs]')
            fig['layout']['yaxis1'].update(title='Cd')
            fig['layout']['yaxis2'].update(title='Cr')
            fig['layout']['yaxis3'].update(title='Cd')
            fig['layout']['yaxis4'].update(title='Cr')
            fig['layout']['yaxis5'].update(title='FSCT')
            fig['layout']['yaxis6'].update(title='FSCT')
            # Execute
            plotly.offline.plot(fig, filename='Plots/estparams.html', auto_open=True)
            # plotly.offline.plot(fig, image = 'png', image_filename='smoothermeas', output_type='file', image_width=1920, image_height=1080 , validate=False)

        ########################################################################################
        
        # if plotRef:
        #     # Plot Reference Data
        #     title1 = 'Reference Position [m]'
        #     trace1 = go.Scatter(x=refMJD, y=refposvel[:,0], mode = 'markers', name = 'x', marker = dict(color = 'rgb(51, 153, 255)'))
        #     trace2 = go.Scatter(x=refMJD, y=refposvel[:,1], mode = 'markers', name = 'y', marker = dict(color = 'rgb(255, 192, 56)'))
        #     trace3 = go.Scatter(x=refMJD, y=refposvel[:,2], mode = 'markers', name = 'z', marker = dict(color = 'rgb(0, 204, 0)'))
        #     trace4 = go.Scatter(x=obsMJD[0,:], y=refobs[:,0], mode = 'markers', name = 'x interp', marker = dict(color = 'rgb(0, 0, 0)'))
        #     trace5 = go.Scatter(x=obsMJD[0,:], y=refobs[:,1], mode = 'markers', name = 'y interp', marker = dict(color = 'rgb(0, 0, 0)'))
        #     trace6 = go.Scatter(x=obsMJD[0,:], y=refobs[:,2], mode = 'markers', name = 'z interp', marker = dict(color = 'rgb(0, 0, 0)'))
        #     title2 = 'Reference Velocity [m/s]'
        #     trace7 = go.Scatter(x=refMJD, y=refposvel[:,3], mode = 'markers', name = 'vx', marker = dict(color = 'rgb(51, 153, 255)'))
        #     trace8 = go.Scatter(x=refMJD, y=refposvel[:,4], mode = 'markers', name = 'vy', marker = dict(color = 'rgb(255, 192, 56)'))
        #     trace9 = go.Scatter(x=refMJD, y=refposvel[:,5], mode = 'markers', name = 'vz', marker = dict(color = 'rgb(0, 204, 0)'))
        #     trace10 = go.Scatter(x=obsMJD[0,:], y=refobs[:,3], mode = 'markers', name = 'vx interp', marker = dict(color = 'rgb(0, 0, 0)'))
        #     trace11 = go.Scatter(x=obsMJD[0,:], y=refobs[:,4], mode = 'markers', name = 'vy interp', marker = dict(color = 'rgb(0, 0, 0)'))
        #     trace12 = go.Scatter(x=obsMJD[0,:], y=refobs[:,5], mode = 'markers', name = 'vz interp', marker = dict(color = 'rgb(0, 0, 0)'))
        #     # Create subplot
        #     fig = tools.make_subplots(rows=3, cols=2, subplot_titles=(title1, title2), shared_xaxes=True)
        #     # Assign traces to subplots
        #     fig.append_trace(trace1, 1, 1)
        #     fig.append_trace(trace2, 2, 1)
        #     fig.append_trace(trace3, 3, 1)
        #     fig.append_trace(trace4, 1, 1)
        #     fig.append_trace(trace5, 2, 1)
        #     fig.append_trace(trace6, 3, 1)

        #     fig.append_trace(trace7, 1, 2)
        #     fig.append_trace(trace8, 2, 2)
        #     fig.append_trace(trace9, 3, 2)
        #     fig.append_trace(trace10, 1, 2)
        #     fig.append_trace(trace11, 2, 2)
        #     fig.append_trace(trace12, 3, 2)
            
        #     # Configure Layout
        #     fig['layout'].update(title='Reference Satellite Data')
        #     # Execute
        #     plotly.offline.plot(fig, filename='reference.html', auto_open=True)

        #####################################################################################################

        TITLE = 'Reference - Filter Position and Velocity Residuals'
        # Plot Reference - Filter Differences
        title1 = 'Position x: Mean = ' + format(meanstdPVfilter[0],'.2f') + ' | STD = ' + format(meanstdPVfilter[1],'.2f')
        trace1 = go.Scatter(x=tstr, y=PVDiff[:,0], mode = 'markers', name = 'x', marker = dict(color = 'rgb(51, 153, 255)'))
        title2 = 'Position y: Mean = ' + format(meanstdPVfilter[2],'.2f') + ' | STD = ' + format(meanstdPVfilter[3],'.2f')
        trace2 = go.Scatter(x=tstr, y=PVDiff[:,1], mode = 'markers', name = 'y', marker = dict(color = 'rgb(255, 192, 56)'))
        title3 = 'Position z: Mean = ' + format(meanstdPVfilter[4],'.2f') + ' | STD = ' + format(meanstdPVfilter[5],'.2f')
        trace3 = go.Scatter(x=tstr, y=PVDiff[:,2], mode = 'markers', name = 'z', marker = dict(color = 'rgb(0, 204, 0)'))
        title4 = 'Velocity vx: Mean = ' + format(meanstdPVfilter[6],'.2f') + ' | STD = ' + format(meanstdPVfilter[7],'.2f')
        trace4 = go.Scatter(x=tstr, y=PVDiff[:,3], mode = 'markers', name = 'vx', marker = dict(color = 'rgb(51, 153, 255)'))
        title5 = 'Velocity vy: Mean = ' + format(meanstdPVfilter[8],'.2f') + ' | STD = ' + format(meanstdPVfilter[9],'.2f')
        trace5 = go.Scatter(x=tstr, y=PVDiff[:,4], mode = 'markers', name = 'vy', marker = dict(color = 'rgb(255, 192, 56)'))
        title6 = 'Velocity vz: Mean = ' + format(meanstdPVfilter[10],'.2f') + ' | STD = ' + format(meanstdPVfilter[11],'.2f')
        trace6 = go.Scatter(x=tstr, y=PVDiff[:,5], mode = 'markers', name = 'vz', marker = dict(color = 'rgb(0, 204, 0)'))
        # Plot 3-sigma Bounds
        trace7 = go.Scatter(x=tstr, y=sigPf[0,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds')
        trace8 = go.Scatter(x=tstr, y=sigPf[1,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace9 = go.Scatter(x=tstr, y=sigPf[2,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace10 = go.Scatter(x=tstr, y=-sigPf[0,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace11 = go.Scatter(x=tstr, y=-sigPf[1,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace12 = go.Scatter(x=tstr, y=-sigPf[2,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace13 = go.Scatter(x=tstr, y=sigPf[3,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend=False)
        trace14 = go.Scatter(x=tstr, y=sigPf[4,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace15 = go.Scatter(x=tstr, y=sigPf[5,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace16 = go.Scatter(x=tstr, y=-sigPf[3,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace17 = go.Scatter(x=tstr, y=-sigPf[4,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace18 = go.Scatter(x=tstr, y=-sigPf[5,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        # Create subplot
        fig = make_subplots(rows=3, cols=2, subplot_titles=(title1, title4, title2, title5, title3, title6), shared_xaxes=True)
        # Assign traces to subplots
        fig.append_trace(trace1, 1, 1)
        fig.append_trace(trace2, 2, 1)
        fig.append_trace(trace3, 3, 1)
        fig.append_trace(trace4, 1, 2)
        fig.append_trace(trace5, 2, 2)
        fig.append_trace(trace6, 3, 2)
        fig.append_trace(trace7, 1, 1)
        fig.append_trace(trace8, 2, 1)
        fig.append_trace(trace9, 3, 1)
        fig.append_trace(trace10, 1, 1)
        fig.append_trace(trace11, 2, 1)
        fig.append_trace(trace12, 3, 1)
        fig.append_trace(trace13, 1, 2)
        fig.append_trace(trace14, 2, 2)
        fig.append_trace(trace15, 3, 2)
        fig.append_trace(trace16, 1, 2)
        fig.append_trace(trace17, 2, 2)
        fig.append_trace(trace18, 3, 2)
        # Configure Layout
        fig['layout'].update(title=TITLE)
        fig['layout']['xaxis1'].update(title='Time since epoch [hrs]')
        fig['layout']['xaxis2'].update(title='Time since epoch [hrs]')
        fig['layout']['yaxis1'].update(title='x [m]')
        fig['layout']['yaxis2'].update(title='vx [m/s]')
        fig['layout']['yaxis3'].update(title='y [m]')
        fig['layout']['yaxis4'].update(title='vy [m/s]')
        fig['layout']['yaxis5'].update(title='z [m]')
        fig['layout']['yaxis6'].update(title='vz [m/s]')
        # Execute
        plotly.offline.plot(fig, filename='Plots/filterresiduals.html', auto_open=True)
        # plotly.offline.plot(fig, image = 'png', image_filename='filterresiduals', output_type='file', image_width=1920, image_height=1080 , validate=False)

        #####################################################################################################

        if smootherRun:
            TITLE = 'Reference - Smoother Position and Velocity Residuals'
            # Plot Reference - Smoother Differences
            title1 = 'Position x: Mean = ' + format(meanstdPVsmooth[0],'.2f') + ' | STD = ' + format(meanstdPVsmooth[1],'.2f')
            trace1 = go.Scatter(x=tstr, y=posvelDiff[:,0], mode = 'markers', name = 'x', marker = dict(color = 'rgb(51, 153, 255)'))
            title2 = 'Position y: Mean = ' + format(meanstdPVsmooth[2],'.2f') + ' | STD = ' + format(meanstdPVsmooth[3],'.2f')
            trace2 = go.Scatter(x=tstr, y=posvelDiff[:,1], mode = 'markers', name = 'y', marker = dict(color = 'rgb(255, 192, 56)'))
            title3 = 'Position z: Mean = ' + format(meanstdPVsmooth[4],'.2f') + ' | STD = ' + format(meanstdPVsmooth[5],'.2f')
            trace3 = go.Scatter(x=tstr, y=posvelDiff[:,2], mode = 'markers', name = 'z', marker = dict(color = 'rgb(0, 204, 0)'))
            title4 = 'Velocity vx: Mean = ' + format(meanstdPVsmooth[6],'.2f') + ' | STD = ' + format(meanstdPVsmooth[7],'.2f')
            trace4 = go.Scatter(x=tstr, y=posvelDiff[:,3], mode = 'markers', name = 'vx', marker = dict(color = 'rgb(51, 153, 255)'))
            title5 = 'Velocity vy: Mean = ' + format(meanstdPVsmooth[8],'.2f') + ' | STD = ' + format(meanstdPVsmooth[9],'.2f')
            trace5 = go.Scatter(x=tstr, y=posvelDiff[:,4], mode = 'markers', name = 'vy', marker = dict(color = 'rgb(255, 192, 56)'))
            title6 = 'Velocity vz: Mean = ' + format(meanstdPVsmooth[10],'.2f') + ' | STD = ' + format(meanstdPVsmooth[11],'.2f')
            trace6 = go.Scatter(x=tstr, y=posvelDiff[:,5], mode = 'markers', name = 'vz', marker = dict(color = 'rgb(0, 204, 0)'))
            # Plot 3-sigma Bounds
            trace7 = go.Scatter(x=tstr, y=sigPs[0,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds')
            trace8 = go.Scatter(x=tstr, y=sigPs[1,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace9 = go.Scatter(x=tstr, y=sigPs[2,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace10 = go.Scatter(x=tstr, y=-sigPs[0,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace11 = go.Scatter(x=tstr, y=-sigPs[1,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace12 = go.Scatter(x=tstr, y=-sigPs[2,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace13 = go.Scatter(x=tstr, y=sigPs[3,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend=False)
            trace14 = go.Scatter(x=tstr, y=sigPs[4,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace15 = go.Scatter(x=tstr, y=sigPs[5,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace16 = go.Scatter(x=tstr, y=-sigPs[3,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace17 = go.Scatter(x=tstr, y=-sigPs[4,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace18 = go.Scatter(x=tstr, y=-sigPs[5,:], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            # Create subplot
            fig = make_subplots(rows=3, cols=2, subplot_titles=(title1, title4, title2, title5, title3, title6), shared_xaxes=True)
            # Assign traces to subplots
            fig.append_trace(trace1, 1, 1)
            fig.append_trace(trace2, 2, 1)
            fig.append_trace(trace3, 3, 1)
            fig.append_trace(trace4, 1, 2)
            fig.append_trace(trace5, 2, 2)
            fig.append_trace(trace6, 3, 2)
            fig.append_trace(trace7, 1, 1)
            fig.append_trace(trace8, 2, 1)
            fig.append_trace(trace9, 3, 1)
            fig.append_trace(trace10, 1, 1)
            fig.append_trace(trace11, 2, 1)
            fig.append_trace(trace12, 3, 1)
            fig.append_trace(trace13, 1, 2)
            fig.append_trace(trace14, 2, 2)
            fig.append_trace(trace15, 3, 2)
            fig.append_trace(trace16, 1, 2)
            fig.append_trace(trace17, 2, 2)
            fig.append_trace(trace18, 3, 2)
            # Configure Layout
            fig['layout'].update(title= TITLE)
            fig['layout']['xaxis1'].update(title='Time since epoch [hrs]')
            fig['layout']['xaxis2'].update(title='Time since epoch [hrs]')
            fig['layout']['yaxis1'].update(title='x [m]')
            fig['layout']['yaxis2'].update(title='vx [m/s]')
            fig['layout']['yaxis3'].update(title='y [m]')
            fig['layout']['yaxis4'].update(title='vy [m/s]')
            fig['layout']['yaxis5'].update(title='z [m]')
            fig['layout']['yaxis6'].update(title='vz [m/s]')
            # Execute
            plotly.offline.plot(fig, filename='Plots/smootherresiduals.html', auto_open=True)
            # plotly.offline.plot(fig, image = 'png', image_filename='smootherresiduals', output_type='file', image_width=1920, image_height=1080 , validate=False)

        #####################################################################################################

        TITLE = 'RIC Reference - Filter Position and Velocity Residuals'
        # Plot RIC Filtered Positions
        title1 = 'Position Radial: Mean = ' + format(meanstdRICfilter[0],'.2f') + ' | STD = ' + format(meanstdRICfilter[1],'.2f')
        title2 = 'Position In-track: Mean = ' + format(meanstdRICfilter[2],'.2f') + ' | STD = ' + format(meanstdRICfilter[3],'.2f')
        title3 = 'Position Cross-track: Mean = ' + format(meanstdRICfilter[4],'.2f') + ' | STD = ' + format(meanstdRICfilter[5],'.2f')
        title4 = 'Velocity Radial: Mean = ' + format(meanstdRICfilter[6],'.2f') + ' | STD = ' + format(meanstdRICfilter[7],'.2f')
        title5 = 'Velocity In-track: Mean = ' + format(meanstdRICfilter[8],'.2f') + ' | STD = ' + format(meanstdRICfilter[9],'.2f')
        title6 = 'Velocity Cross-track: Mean = ' + format(meanstdRICfilter[10],'.2f') + ' | STD = ' + format(meanstdRICfilter[11],'.2f')
        trace1 = go.Scatter(x=tstr, y=xDiffRIC[:,0], mode = 'markers', name = 'radial', marker = dict(color = 'rgb(51, 153, 255)'))
        trace2 = go.Scatter(x=tstr, y=xDiffRIC[:,1], mode = 'markers', name = 'in-track', marker = dict(color = 'rgb(255, 192, 56)'))
        trace3 = go.Scatter(x=tstr, y=xDiffRIC[:,2], mode = 'markers', name = 'cross-track', marker = dict(color = 'rgb(0, 204, 0)'))
        trace4 = go.Scatter(x=tstr, y=xDiffRIC[:,3], mode = 'markers', name = 'radial', marker = dict(color = 'rgb(51, 153, 255)'), showlegend = False)
        trace5 = go.Scatter(x=tstr, y=xDiffRIC[:,4], mode = 'markers', name = 'in-track', marker = dict(color = 'rgb(255, 192, 56)'), showlegend = False)
        trace6 = go.Scatter(x=tstr, y=xDiffRIC[:,5], mode = 'markers', name = 'cross-track', marker = dict(color = 'rgb(0, 204, 0)'), showlegend = False)
        # Plot 3-sigma Bounds
        trace7 = go.Scatter(x=tstr, y=sigPRIC[0], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds')
        trace8 = go.Scatter(x=tstr, y=sigPRIC[1], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace9 = go.Scatter(x=tstr, y=sigPRIC[2], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace10 = go.Scatter(x=tstr, y=-sigPRIC[0], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace11 = go.Scatter(x=tstr, y=-sigPRIC[1], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace12 = go.Scatter(x=tstr, y=-sigPRIC[2], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace13 = go.Scatter(x=tstr, y=sigPRIC[3], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace14 = go.Scatter(x=tstr, y=sigPRIC[4], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace15 = go.Scatter(x=tstr, y=sigPRIC[5], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace16 = go.Scatter(x=tstr, y=-sigPRIC[3], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace17 = go.Scatter(x=tstr, y=-sigPRIC[4], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        trace18 = go.Scatter(x=tstr, y=-sigPRIC[5], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
        # Create subplot
        fig = make_subplots(rows=3, cols=2, subplot_titles=(title1, title4, title2, title5, title3, title6), shared_xaxes=True)
        # Assign traces to subplots
        fig.append_trace(trace1, 1, 1)
        fig.append_trace(trace2, 2, 1)
        fig.append_trace(trace3, 3, 1)
        fig.append_trace(trace4, 1, 2)
        fig.append_trace(trace5, 2, 2)
        fig.append_trace(trace6, 3, 2)
        fig.append_trace(trace7, 1, 1)
        fig.append_trace(trace8, 2, 1)
        fig.append_trace(trace9, 3, 1)
        fig.append_trace(trace10, 1, 1)
        fig.append_trace(trace11, 2, 1)
        fig.append_trace(trace12, 3, 1)
        fig.append_trace(trace13, 1, 2)
        fig.append_trace(trace14, 2, 2)
        fig.append_trace(trace15, 3, 2)
        fig.append_trace(trace16, 1, 2)
        fig.append_trace(trace17, 2, 2)
        fig.append_trace(trace18, 3, 2)
        # Configure Layout
        fig['layout'].update(title= TITLE)
        fig['layout']['xaxis1'].update(title='Time since epoch [hrs]')
        fig['layout']['xaxis2'].update(title='Time since epoch [hrs]')
        fig['layout']['yaxis1'].update(title='Radial [m]')
        fig['layout']['yaxis2'].update(title='Radial [m/s]')
        fig['layout']['yaxis3'].update(title='In-track [m]')
        fig['layout']['yaxis4'].update(title='In-track [m/s]')
        fig['layout']['yaxis5'].update(title='Cross-track [m]')
        fig['layout']['yaxis6'].update(title='Cross-track [m/s]')
        # Execute
        plotly.offline.plot(fig, filename='Plots/RICfilterresiduals.html', auto_open=True)
        # plotly.offline.plot(fig, image = 'png', image_filename='RICfilterresiduals', output_type='file', image_width=1920, image_height=1080 , validate=False)

        #####################################################################################################

        if smootherRun:
            TITLE = 'RIC Reference - Smoother Position and Velocity Residuals'
            # Plot RIC Smoother Positions
            title1 = 'Position Radial: Mean = ' + format(meanstdRICsmooth[0],'.2f') + ' | STD = ' + format(meanstdRICsmooth[1],'.2f')
            title2 = 'Position In-track: Mean = ' + format(meanstdRICsmooth[2],'.2f') + ' | STD = ' + format(meanstdRICsmooth[3],'.2f')
            title3 = 'Position Cross-track: Mean = ' + format(meanstdRICsmooth[4],'.2f') + ' | STD = ' + format(meanstdRICsmooth[5],'.2f')
            title4 = 'Velocity Radial: Mean = ' + format(meanstdRICsmooth[6],'.2f') + ' | STD = ' + format(meanstdRICsmooth[7],'.2f')
            title5 = 'Velocity In-track: Mean = ' + format(meanstdRICsmooth[8],'.2f') + ' | STD = ' + format(meanstdRICsmooth[9],'.2f')
            title6 = 'Velocity Cross-track: Mean = ' + format(meanstdRICsmooth[10],'.2f') + ' | STD = ' + format(meanstdRICsmooth[11],'.2f')
            trace1 = go.Scatter(x=tstr, y=xSDiffRIC[:,0], mode = 'markers', name = 'radial', marker = dict(color = 'rgb(51, 153, 255)'))
            trace2 = go.Scatter(x=tstr, y=xSDiffRIC[:,1], mode = 'markers', name = 'in-track', marker = dict(color = 'rgb(255, 192, 56)'))
            trace3 = go.Scatter(x=tstr, y=xSDiffRIC[:,2], mode = 'markers', name = 'cross-track', marker = dict(color = 'rgb(0, 204, 0)'))
            trace4 = go.Scatter(x=tstr, y=xSDiffRIC[:,3], mode = 'markers', name = 'radial', marker = dict(color = 'rgb(51, 153, 255)'), showlegend = False)
            trace5 = go.Scatter(x=tstr, y=xSDiffRIC[:,4], mode = 'markers', name = 'in-track', marker = dict(color = 'rgb(255, 192, 56)'), showlegend = False)
            trace6 = go.Scatter(x=tstr, y=xSDiffRIC[:,5], mode = 'markers', name = 'cross-track', marker = dict(color = 'rgb(0, 204, 0)'), showlegend = False)
            # Plot 3-sigma Bounds
            trace7 = go.Scatter(x=tstr, y=sigPSRIC[0], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds')
            trace8 = go.Scatter(x=tstr, y=sigPSRIC[1], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace9 = go.Scatter(x=tstr, y=sigPSRIC[2], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace10 = go.Scatter(x=tstr, y=-sigPSRIC[0], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace11 = go.Scatter(x=tstr, y=-sigPSRIC[1], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace12 = go.Scatter(x=tstr, y=-sigPSRIC[2], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace13 = go.Scatter(x=tstr, y=sigPSRIC[3], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace14 = go.Scatter(x=tstr, y=sigPSRIC[4], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace15 = go.Scatter(x=tstr, y=sigPSRIC[5], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace16 = go.Scatter(x=tstr, y=-sigPSRIC[3], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace17 = go.Scatter(x=tstr, y=-sigPSRIC[4], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            trace18 = go.Scatter(x=tstr, y=-sigPSRIC[5], line = dict(color = ('rgb(255, 51, 51)')), name = '3-sigma bounds', showlegend = False)
            # Create subplot
            fig = make_subplots(rows=3, cols=2, subplot_titles=(title1, title4, title2, title5, title3, title6), shared_xaxes=True)
            # Assign traces to subplots
            fig.append_trace(trace1, 1, 1)
            fig.append_trace(trace2, 2, 1)
            fig.append_trace(trace3, 3, 1)
            fig.append_trace(trace4, 1, 2)
            fig.append_trace(trace5, 2, 2)
            fig.append_trace(trace6, 3, 2)
            fig.append_trace(trace7, 1, 1)
            fig.append_trace(trace8, 2, 1)
            fig.append_trace(trace9, 3, 1)
            fig.append_trace(trace10, 1, 1)
            fig.append_trace(trace11, 2, 1)
            fig.append_trace(trace12, 3, 1)
            fig.append_trace(trace13, 1, 2)
            fig.append_trace(trace14, 2, 2)
            fig.append_trace(trace15, 3, 2)
            fig.append_trace(trace16, 1, 2)
            fig.append_trace(trace17, 2, 2)
            fig.append_trace(trace18, 3, 2)
            # Configure Layout
            fig['layout'].update(title=TITLE)
            fig['layout']['xaxis1'].update(title='Time since epoch [hrs]')
            fig['layout']['xaxis2'].update(title='Time since epoch [hrs]')
            fig['layout']['yaxis1'].update(title='Radial [m]')
            fig['layout']['yaxis2'].update(title='Radial [m/s]')
            fig['layout']['yaxis3'].update(title='In-track [m]')
            fig['layout']['yaxis4'].update(title='In-track [m/s]')
            fig['layout']['yaxis5'].update(title='Cross-track [m]')
            fig['layout']['yaxis6'].update(title='Cross-track [m/s]')
            # Execute
            plotly.offline.plot(fig, filename='Plots/RICsmootherresiduals.html', auto_open=True)
            # plotly.offline.plot(fig, image = 'png', image_filename='RICsmootherresiduals', output_type='file', image_width=1920, image_height=1080 , validate=False)

        #####################################################################################################