from datetime import datetime, timedelta
import numpy as np
import scipy.linalg
import scipy.io as sio
from scipy.interpolate import interp1d
import json
import math
import os
import sys
import time
import math
import plotly
import plotly.graph_objs as go
from plotly import tools
from astropy.time import Time
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from orbdetpy import cart2rangerate
from orbdetpy import cart2range
from orbdetpy import cart2azel
from orbdetpy import cart2radec
from orbdetpy import lla2eciJ2000

def matlab_to_python_datetime(matlab_datenum):
    """Convert matlab datenumes into date time"""
    return datetime.fromordinal(int(matlab_datenum)) + timedelta(days=matlab_datenum % 1) - timedelta(days=366)


def sigmapoints(xBar, PBar):
    """Unscaled Unscented Transform with no mean"""
    n = np.size(xBar, 0)
    L = scipy.linalg.cholesky(PBar, lower=True)
    sigmaPoints = np.zeros([n, 2*n])
    sigmaPoints[:, 0:n] = xBar + np.sqrt(n) * L
    sigmaPoints[:, n:2*n] = xBar - np.sqrt(n) * L
    return sigmaPoints


def sigmapoints2(xBar, PBar):
    """Unscaled Unscented Transform (include mean)"""
    n = np.size(xBar, 0)
    L = scipy.linalg.cholesky(PBar, lower=True)
    sigmaPoints = np.zeros([n, 2*n+1])
    sigmaPoints[:, [0]] = xBar
    sigmaPoints[:, 1:n+1] = xBar + np.sqrt(n) * L
    sigmaPoints[:, n+1:2*n+1] = xBar - np.sqrt(n) * L
    return sigmaPoints


def sigmapoints3(xBar, PBar, alpha, k):
    """Scaled Unscented Transform (include mean)"""
    n = np.size(xBar, 0)
    lambda_ = alpha ** 2 * (n + k) - n
    L = scipy.linalg.cholesky(PBar, lower=True)
    sigmaPoints = np.empty([n, 2 * n])
    sigmaPoints[:, [0]] = xBar
    sigmaPoints[:, 1:n + 1] = xBar + np.sqrt(n + lambda_) * L
    sigmaPoints[:, n:2 * n + 1] = xBar - np.sqrt(n + lambda_) * L
    return sigmaPoints


def readODoutJSON(filename):
    """Read in JSON file output from testodet.py"""
    with open(filename) as f:
        jsonData = json.load(f)

    dataSize0 = np.size(jsonData["Estimation"])
    dataSize1 = np.size(jsonData["Estimation"][0]["xBar"])
    xPre = np.zeros([dataSize1, dataSize0])
    xPost = np.zeros([dataSize1, dataSize0])
    PPre = np.zeros([dataSize1, dataSize1, dataSize0])
    PPre2 = []
    PPost = np.zeros([dataSize1, dataSize1, dataSize0])
    PPost2 = []
    SPpre = np.zeros([dataSize1, dataSize1*2, dataSize0])
    SPpre2 = []

    for ii in range(0, dataSize0, 1):
        xPre[:, [ii]] = np.reshape(np.asarray(
            jsonData["Estimation"][ii]["xBar"]), [dataSize1, 1])
        xPost[:, [ii]] = np.reshape(np.asarray(
            jsonData["Estimation"][ii]["EstimatedState"]), [dataSize1, 1])

        PPre2 = np.array(jsonData["Estimation"][ii]["PBar"])
        PPre2 = np.reshape(PPre2.T, [dataSize1, dataSize1, 1])
        PPre[:, :, [ii]] = PPre2

        PPost2 = np.array(jsonData["Estimation"][ii]["EstimatedCovariance"])
        PPost2 = np.reshape(PPost2.T, [dataSize1, dataSize1, 1])
        PPost[:, :, [ii]] = PPost2

        SPpre2 = np.array(jsonData["Estimation"][ii]["SPpre"])
        SPpre2 = np.reshape(SPpre2, [dataSize1, dataSize1*2, 1])
        SPpre[:, :, [ii]] = SPpre2

    return xPre, xPost, PPre, PPost, SPpre


def dataConfigMAT(filename):

    filestore = sio.loadmat(filename)

    xPre = filestore['xfm']
    xPost = filestore['xfp']
    PPre = filestore['Pfm']
    PPost = filestore['Pfp']
    SPpre = filestore['Xfm']

    xS = filestore['xs']
    PS = filestore['Ps']

    RS = filestore['Rs']

    return xPre, xPost, PPre, PPost, SPpre, xS, PS, RS


def interpolate1d(time, newtime, posvel):
    """interpolated using python interp1d"""
    # interpolates a given cartesian [6x1] at time to newtime
    fx = interp1d(time, posvel[:, 0], kind='cubic')
    fy = interp1d(time, posvel[:, 1], kind='cubic')
    fz = interp1d(time, posvel[:, 2], kind='cubic')
    fvx = interp1d(time, posvel[:, 3], kind='cubic')
    fvy = interp1d(time, posvel[:, 4], kind='cubic')
    fvz = interp1d(time, posvel[:, 5], kind='cubic')

    indexsize = np.size(newtime)
    PVobs = np.zeros([indexsize, 6])
    PVobs[:, 0] = fx(newtime)
    PVobs[:, 1] = fy(newtime)
    PVobs[:, 2] = fz(newtime)
    PVobs[:, 3] = fvx(newtime)
    PVobs[:, 4] = fvy(newtime)
    PVobs[:, 5] = fvz(newtime)

    return PVobs


def lla2xyz(lat, lon, alt):
    """WGS84 planet fixed lat, lon, alt [rad, rad, km] to x, y, z [km]"""
    dtr = math.pi / 180
    a = 6378.1370
    b = 6356.7523
    f = (a-b) / a
    lat = lat * dtr
    lon = lon * dtr
    sinlat = math.sin(lat)
    coslat = math.cos(lat)
    sinlat2 = sinlat * sinlat
    denom = np.sqrt(1 - (2*f - f*f)*sinlat2)
    G1 = a / denom + alt
    G2 = a*(1-f)**2 / denom + alt

    x = G1 * coslat * math.cos(lon)
    y = G1 * coslat * math.sin(lon)
    z = G2 * sinlat

    return x, y, z


def xyz2lla(x, y, z):
    """WGS84 planet fixed x, y, z [km] to lat, lon, alt [rad, rad, km]"""
    rtd = 180 / math.pi

    a = 6378.1370
    b = 6356.7523

    a2 = a * a
    b2 = b * b

    f = (a - b) / a
    e2 = f * (2 - f)
    ep2 = (a2 - b2) / b2
    p = np.sqrt(x**2 + y**2)

    # Two sides and hypotenuse of right-angle-triangle with one angle=theta
    s1 = z * a
    s2 = p * b
    h = np.sqrt(s1**2 + s2**2)
    sin_theta = s1 / h
    cos_theta = s2 / h
    theta = math.atan(s1 / s2)

    # Two sides and hypotenuse of right-angle-triangle with one angle=lat
    s1 = z + ep2 * b * (sin_theta**3)
    s2 = p - a * e2 * (cos_theta**3)
    h = math.sqrt(s1**2 + s2**2)

    tan_lat = s1 / s2
    sin_lat = s1 / h
    cos_lat = s2 / h
    lat = math.atan(tan_lat)
    lat = rtd * lat

    #N = a2 * ((a2 * (cos_lat**2) + b2 * (sin_lat**2))**(-0.5))
    N = a2 / math.sqrt(a2 * (cos_lat**2) + b2 * (sin_lat**2))
    alt = p / cos_lat - N

    lon = rtd * math.atan2(y, x)

    return lat, lon, alt


def xECI2RIC(xECI, r, v):
    """Transform ECI vector to RIC
    :param xECI: 3x1 pos [km] 
    :param r: 3x1 reference ECI pos [km]
    :param v: 3x1 reference ECI vel [km/s]
    """
    xECI = np.reshape(xECI, [3, 1])
    r = np.reshape(r, [3, 1])
    v = np.reshape(v, [3, 1])

    nr = np.linalg.norm(r)
    ruv = r / nr
    crossrv = np.cross(r, v, axis=0)
    cuv = crossrv / np.linalg.norm(crossrv)
    auv = np.cross(cuv, ruv, axis=0)
    A = np.zeros([3, 3])
    A[0, :] = ruv.T
    A[1, :] = auv.T
    A[2, :] = cuv.T
    xRIC = np.matmul(A, xECI)
    xRIC = np.squeeze(xRIC)

    return xRIC


def PECI2RIC(PECI, r, v):
    """Transform ECI vector to RIC
    :param PECI: 6x6 pos [km] vel [km/s] 
    :param r: 3x1 reference ECI pos [km]
    :param v: 3x1 reference ECI vel [km/s]
    """
    PECI = np.reshape(PECI, [6, 6])
    r = np.reshape(r, [3, 1])
    v = np.reshape(v, [3, 1])
    nr = np.linalg.norm(r)
    ruv = r / nr
    crossrv = np.cross(r, v, axis=0)
    cuv = crossrv / np.linalg.norm(crossrv)
    crosscr = np.cross(cuv, ruv, axis=0)
    auv = crosscr / np.linalg.norm(crosscr)
    A = np.zeros([3, 3])
    A[0, :] = ruv.T
    A[1, :] = auv.T
    A[2, :] = cuv.T
    C = np.zeros([6, 6])
    C[0:3, 0:3] = A
    C[3:6, 3:6] = A
    PRIC = np.matmul(C, np.matmul(PECI, C.T))

    return PRIC


def writeObsJSON(fnObs, startMeasurement, endMeasurement, obs_isot, PVobs, RA, DEC):
    """Write .json file from matlab data"""
    output = '[\n  '
    for ii in range(startMeasurement, endMeasurement, 1):
        obsTime = obs_isot[ii]+'Z'
        output = output + '{\n    "TrueState": {\n      "Cartesian": [\n        '+str(PVobs[ii, 0])+',\n        '+str(PVobs[ii, 1])+',\n        '+str(
            PVobs[ii, 2])+',\n        '+str(PVobs[ii, 3])+',\n        '+str(PVobs[ii, 4])+',\n        '+str(PVobs[ii, 5])+'\n      ]\n    },\n    '
        output = output + '"Time": "'+obsTime+'",\n    '
        output = output + '"Station": "Station1",\n    '
        output = output + '"RightAscension": ' + \
            str(float(RA[0, ii]*math.pi/180))+',\n    '
        output = output + '"Declination": ' + \
            str(float(DEC[0, ii]*math.pi/180))+'\n  },\n  '

    output = output[:-4] + '\n]'

    with open(fnObs, "w") as fp:
        fp.write(output)

    return


def rtsSmoother(xPre, xPost, PPre, PPost, SPpre):
    """ 
    Runs an Unscented Rauch-Tung-Striebel Smoother
    :param xPre: Array of filtered means at each time step (pre-update) 
    :param xPost: Array of filtered means at each time step (post-update)
    :param PPre: Array of filtered covariancess at each time step (pre-update), includes process noise
    :param PPost: Array of filtered covariancess at each time step (post-update)
    :param SPpre: Array of filtered SPs at each time step (pre-update), Does not include process noise, Unscaled Transform, does not include mean
    :type xPre: np.array [statelength x number of observations] 
    :type xPost: np.array [statelength x number of observations] 
    :type PPre: np.array [statelength x statelength x number of observations]
    :type PPost: np.array [statelength x statelength x number of observations]
    :type SPpre: np.array [statelength x 2*statelength x number of observations]
    :return xS: Array of smoothed means at each time step
    :return PS: Array of smoothed covariances at each time step
    :return Rs: McReynolds Consistency Metric (consistency between filtered and smoothed solution)
    :return time: array from 1 to number of observations
    :rtype xS: np.array [statelength x number of observations] 
    :rtype PS: np.array [statelength x statelength x number of observations]
    :rtype Rs: np.array [statelength x number of observations] 
    :rytpe time: linspace [number of observations]
     """
    # Unscaled UKF RTS Smoother
    N = np.size(xPost, 0)
    T = np.size(xPost, 1)
    time = np.linspace(1, T, T)
    xS = np.zeros([N, T])
    PS = np.zeros([N, N, T])
    Rs = np.zeros([N, T])
    T = T - 1

    wi_c = 1/(2*N)
    xS[:, [T]] = xPost[:, [T]]
    PS[:, :, T] = PPost[:, :, T]

    for k in range(T-1, -1, -1):
        # Extract filter means and covariances
        xkp1 = xPre[:, [k+1]]
        Pkp1 = PPre[:, :, k+1]
        Xkp1 = SPpre[:, :, k+1]
        xk = xPost[:, [k]]
        Pk = PPost[:, :, k]

        # Compute sigmapoints from state and covariance, (take in pre filtered sigmapoints directly)
        #SPpre = sigmapoints2(xPre[:, [kk]], PPre[:, :, kk])
        Xk = sigmapoints(xk, Pk)

        # Compute cross-covariance between t = k and t = k+1
        Ckp1 = np.zeros(N)
        for i in range(0, 2*N, 1):
            Ckp1 = Ckp1 + wi_c * (Xk[:, [i]] - xk) * \
                np.transpose(Xkp1[:, [i]] - xkp1)

        # Smoother update
        Ak = np.dot(Ckp1, np.linalg.inv(Pkp1))
        xS[:, [k]] = xk + np.dot(Ak, (xS[:, [k+1]] - xkp1))
        PS[:, :, k] = Pk + \
            np.dot(Ak, np.dot((PS[:, :, k+1] - Pkp1), np.transpose(Ak)))

        # McReynold's Filter Smoother Consistency Test
        dx = xk - xS[:, [k]]
        dP = Pk - PS[:, :, k]
        sig = np.sqrt(np.diagonal(dP))
        sig = np.reshape(sig, [N, 1])
        Rs[:, [k]] = np.divide(dx, sig)
        # for j in range(0,N,1):
        #     if abs(Rs[j,[k]]) > 3:
        #         print('Warning, consistency test on obs', k, 'row', j, Rs[j,[k]])

    return xS, PS, Rs, time


def get3sig(cov):
    """ Takes 3 times the square root of the diagonal of a covariance matrix
    :param cov: covariance matrix 
    :type meas: np.array [6x6]
    :return sig: 3*sigma values for each diagonal element
    :type sig: np.array [6x1]
    """
    # cov = [6,6,numObs] pos vel diagonal covariance
    num = cov.shape[2]
    sig = np.zeros([6, num])
    sig[0, :] = 3*np.sqrt(cov[0, 0, :])
    sig[1, :] = 3*np.sqrt(cov[1, 1, :])
    sig[2, :] = 3*np.sqrt(cov[2, 2, :])
    sig[3, :] = 3*np.sqrt(cov[3, 3, :])
    sig[4, :] = 3*np.sqrt(cov[4, 4, :])
    sig[5, :] = 3*np.sqrt(cov[5, 5, :])

    return sig


def writeODcfgJSON(fnODin, fnODout, startTime, endTime, startMeasurement, PVobs, lat, lon, alt):
    """Updates fnODin to the given inputs
    :param startTime: propagation start time
    :param endTime: propagation end time
    :param startMeasurement: measurement number to set initial state
    :param PVobs: position and velocity at all times (only need start measurement)
    :param lat: geodetic latitude of station [rad]
    :param lon: geodetic longitude of station [rad]
    :param alt: altitude of station [m]
    """
    ### OD CONFIG FILE ###
    with open(fnODin) as f:
        jsonData = json.load(f)

        # Configured from input data
        jsonData['Propagation']['Start'] = startTime+'Z'
        jsonData['Propagation']['End'] = endTime+'Z'
        jsonData['Propagation']['InitialState'] = PVobs[startMeasurement, :].tolist()

        jsonData['Stations']['Station1']['Latitude'] = lat
        jsonData['Stations']['Station1']['Longitude'] = lon
        jsonData['Stations']['Station1']['Altitude'] = alt

        # Manually Configure
        jsonData['Drag']['Coefficient']['Estimation'] = 'Estimate'
        jsonData['RadiationPressure']['Creflection']['Estimation'] = 'Estimate'

        jsonData['Measurements']['RightAscension']['Error'] = [
            10 * math.pi / (180*3600)]
        jsonData['Measurements']['Declination']['Error'] = [
            10 * math.pi / (180*3600)]

        sigx = 100 * (1000) ** 2
        sigy = 100 * (1000) ** 2
        sigz = 100 * (1000) ** 2
        sigvx = 0.5 ** 2
        sigvy = 0.5 ** 2
        sigvz = 0.5 ** 2
        jsonData['Estimation']['Covariance'] = [
            sigx, sigy, sigz, sigvx, sigvy, sigvz, 1.00, 0.25]

        sigx = 1E-24
        sigy = 1E-24
        sigz = 1E-24
        sigvx = 1E-24
        sigvy = 1E-24
        sigvz = 1E-24
        jsonData['Estimation']['ProcessNoise'] = [
            sigx, sigy, sigz, sigvx, sigvy, sigvz]

    with open(fnODout, "w") as f:
        json.dump(jsonData, f, indent=4)

    return


def meanstdANG(meas):
    """Gets mean and std from meas
    :param meas: angle1 and angle2 [rad] 
    :type meas: np.array [numObsx2]
    :return meanstdANG: mean and std of meas [4]
    :type meanstdANG: list
    """
    # meas: [numObsx2] RA [rad] DEC [rad]
    meanstdANG = [np.mean(meas[:, 0]), np.std(
        meas[:, 0]), np.mean(meas[:, 1]), np.std(meas[:, 1])]
    # meanstd = [meanRA stdRA meanDEC stdDEC]

    return meanstdANG


def pv2radec(cfg, gslat, gslon, gsalt, time, angular, sigma, pv):
    """Converts a position and velocity into a right ascension & declination measurement
    :param cfg: od cfg filename
    :param gslat: geodetic latitude of station [rad]
    :param gslon: geodetic longitude of station [rad]
    :param gsalt: altitude of station [m]
    :param time: time stamp of observation
    :param angular: observed RA and DEC at time [rad]
    :param sigma: measurement noise of observed RA and DEC [rad]
    :param pv: position and velocity [m]
    :type cfg: str
    :type gslat: float
    :type gslon: float
    :type gsalt: float
    :type time: str
    :type angular: list
    :type sigma: list
    :type pv: list
    :return radec: RA and DEC measurement at time [rad]
    :rtype: list
    """
    with open(cfg, "r") as fp:
        cfgjson = fp.read()

    radec = cart2radec(cfgjson, gslat, gslon, gsalt, time, angular, sigma, pv)

    return(radec)


def pv2azel(cfg, gslat, gslon, gsalt, time, angular, sigma, pv):
    """Converts a position and velocity into a right ascension & declination measurement
    :param cfg: od cfg filename
    :param gslat: geodetic latitude of station [rad]
    :param gslon: geodetic longitude of station [rad]
    :param gsalt: altitude of station [m]
    :param time: time stamp of observation
    :param angular: observed AZ and EL at time [rad]
    :param sigma: measurement noise of observed RA and DEC [rad]
    :param pv: position and velocity [m]
    :type cfg: str
    :type gslat: float
    :type gslon: float
    :type gsalt: float
    :type time: str
    :type angular: list
    :type sigma: list
    :type pv: list
    :return azel: AZ and EL measurement at time [rad]
    :rtype: list
    """
    with open(cfg, "r") as fp:
        cfgjson = fp.read()

    azel = cart2azel(cfgjson, gslat, gslon, gsalt, time, angular, sigma, pv)

    return(azel)


def pv2range(cfg, gslat, gslon, gsalt, time, range1, error, pv):
    """Converts a position and velocity into a right ascension & declination measurement
    :param cfg: od cfg filename
    :param gslat: geodetic latitude of station [rad]
    :param gslon: geodetic longitude of station [rad]
    :param gsalt: altitude of station [m]
    :param time: time stamp of observation
    :param range1: observed range at time [m]
    :param sigma: measurement noise of observed RA and DEC [rad]
    :param pv: position and velocity [m]
    :type cfg: str
    :type gslat: float
    :type gslon: float
    :type gsalt: float
    :type time: str
    :type range1: list
    :type error: list
    :type pv: list
    :return azel: AZ and EL measurement at time [rad]
    :rtype: list
    """
    with open(cfg, "r") as fp:
        cfgjson = fp.read()

    range1 = cart2range(cfgjson, gslat, gslon, gsalt, time, range1, error, pv)

    return(range1)


def pv2rangerate(cfg, gslat, gslon, gsalt, time, rangerate1, error, pv):
    """Converts a position and velocity into a right ascension & declination measurement
    :param cfg: od cfg filename
    :param gslat: geodetic latitude of station [rad]
    :param gslon: geodetic longitude of station [rad]
    :param gsalt: altitude of station [m]
    :param time: time stamp of observation
    :param rangerate1: observed rangerate at time [m]
    :param sigma: measurement noise of observed RA and DEC [rad]
    :param pv: position and velocity [m]
    :type cfg: str
    :type gslat: float
    :type gslon: float
    :type gsalt: float
    :type time: str
    :type range1: list
    :type error: list
    :type pv: list
    :return rangerate: rangerate measurement at time [rad]
    :rtype: list
    """
    with open(cfg, "r") as fp:
        cfgjson = fp.read()

    rangerate1 = cart2rangerate(
        cfgjson, gslat, gslon, gsalt, time, rangerate1, error, pv)

    return(rangerate1)


def meanstdPV(PV):
    """Gets mean and std from PV
    :param PV: position and velocity 
    :type PV: np.array [numObsx6]
    :return meanstdPV: mean and std of PV [12]
    :type meanstdPV: list
    """
    meanstdPV = [np.mean(PV[:, 0]), np.std(PV[:, 0]), np.mean(PV[:, 1]), np.std(PV[:, 1]), np.mean(PV[:, 2]), np.std(
        PV[:, 2]), np.mean(PV[:, 3]), np.std(PV[:, 3]), np.mean(PV[:, 4]), np.std(PV[:, 4]), np.mean(PV[:, 5]), np.std(PV[:, 5])]

    return meanstdPV


def lagrange(x, i, xm):
    """Evaluates the i-th Lagrange polynomial at x based on grid data xm"""
    n = len(xm) - 1
    y = 1.
    for j in range(n+1):
        if i != j:
            y *= (x-xm[j])/(xm[i]-xm[j])
    return y


def interpolation(x, xm, ym):
    """Interpolated grid (xm,ym) points to (x,y)"""
    n = len(xm) - 1
    lagrpoly = np.array([lagrange(x, i, xm) for i in range(n+1)])
    y = np.dot(ym, lagrpoly)
    return y


def interpolationPV(newtime, time, PV, degree):
    """Interpolates (time,PV) to (newtime,newPV)
    :param newtime: desired time to interpolate PV to
    :param time: time associated with PV
    :param PV: position and veloicty at given time
    :param degree: desired lagrange polynomial order 
    """
    indexsize = np.size(time)
    val = int(indexsize/10)
    x = interpolation(newtime, time[0::val], PV[0::val, 0])
    y = interpolation(newtime, time[0::val], PV[0::val, 1])
    z = interpolation(newtime, time[0::val], PV[0::val, 2])
    vx = interpolation(newtime, time[0::val], PV[0::val, 3])
    vy = interpolation(newtime, time[0::val], PV[0::val, 4])
    vz = interpolation(newtime, time[0::val], PV[0::val, 5])

    indexsize = np.size(newtime)
    PVobs = np.zeros([indexsize, 6])
    PVobs[:, 0] = x
    PVobs[:, 1] = y
    PVobs[:, 2] = z
    PVobs[:, 3] = vx
    PVobs[:, 4] = vy
    PVobs[:, 5] = vz

    return PVobs


def interpolationANG(newtime, time, ang, degree):
    """Interpolates (time,PV) to (newtime,newPV)
    :param newtime: desired time to interpolate PV to
    :param time: time associated with PV
    :param PV: position and veloicty at given time
    :param degree: desired lagrange polynomial order 
    """
    indexsize = np.size(time)
    val = int(indexsize/10)
    ang1 = interpolation(newtime, time[0::val], ang[0::val, 0])
    ang2 = interpolation(newtime, time[0::val], ang[0::val, 1])

    indexsize = np.size(newtime)
    ANGobs = np.zeros([indexsize, 2])
    ANGobs[:, 0] = ang1
    ANGobs[:, 1] = ang2

    return ANGobs


def aberrationCorrection(obs_data, plot):
    """Calculate and apply effect of aberration to provided observation file
    :param obs_data: filename of obs data file
    :param plot: 1 or 0 to plot correction values over time
    :type obs_data: str
    :type plot: 1 or 0
    """
    # Reference: pg 143-152 Astronomical Algorithms 2nd ed. Jean Meeus

    with open(obs_data, "r") as fp:
            inp = json.load(fp)

    tstr, ang1, ang2 = [], [], []
    for i in inp:
            # Get time string
            tstr.append(i["Time"])
            # Get obs
            ang1.append(i['RightAscension'])
            ang2.append(i['Declination'])

    time_UTC = Time(tstr, format='isot', scale='utc')
    time_TT = time_UTC.tt
    # JDE = Julian Day Ephemeris
    JDE = time_TT.jd

    # Degrees to radians
    D2R = math.pi/180
    # Radians to degrees
    R2D = 180/math.pi
    # Radians to arcseconds
    R2A = 648000.0/math.pi
    # Arcseconds to radians
    A2R = math.pi/648000.

    numObs = len(tstr)
    devRA, devDEC = [], []
    for ii in range(0, numObs, 1):
        # T = time in Julian Centuries pg 59
        T = (JDE[ii] - 2451545)/36525
        T2 = T**2
        T3 = T**3
        # kappa = constant of aberration [arcsec] pg 151
        kappa = A2R*20.49552
        # L0 = geometric mean longitude of the Sun (mean equinox of the date) [rad] pg 163
        L0 = D2R*(280.46646 + 36000.76983*T + 0.0003032*T2)
        # M = Mean anomaly of the Sun [rad] pg 163
        M = D2R*(357.52911 + 35999.05029*T - 0.0001537*T2)
        # C = Sun's equation of the center [rad] pg 164
        C = D2R*((1.914602-0.004817*T-0.000014*T2)*math.sin(M) +
                 (0.019993-0.000101*T)*math.sin(2*M) + 0.000289*math.sin(3*M))
        # lonS = true (geometric) longitude of the Sun [rad] pg 163
        lonS = L0 + C
        while lonS < 0:
            lonS = lonS + 2*math.pi
        while lonS > 2*math.pi:
            lonS = lonS - 2*math.pi
        # e = eccentricity of the Earth's orbit
        e = 0.016708634 - 0.000042037*T - 0.0000001267*T2
        # pi = longitude of the perihelion [rad]
        pi = D2R*(102.93735 + 1.71946*T + 0.00046*T2)
        # epsilon = obliquity of the ecliptic
        eps = 21.448 - 46.8150*T - 0.00059*T2 + 0.001813*T3
        eps = D2R*(23 + 26/60 + eps/3600)
        # Right Ascension observation [rad]
        alpha = ang1[ii]
        # Declination observation [rad]
        dec = ang2[ii]
        # abbreviate calculations
        CA = math.cos(alpha)
        CS = math.cos(lonS)
        CD = math.cos(dec)
        CP = math.cos(pi)
        CE = math.cos(eps)
        SA = math.sin(alpha)
        SS = math.sin(lonS)
        SP = math.sin(pi)
        SD = math.sin(dec)
        TE = math.tan(eps)
        # dalpha = deviation in right ascension [rad]
        dalpha = -kappa*(CA*CS*CE + SA*SS)/CD + e*kappa*(CA*CP*CE + SA*SP)/CD
        # ddec = deviation in declination [rad]
        ddec = -kappa*(CS*CE*(TE*CD - SA*SD)+CA*SD*SS) + e * \
            kappa*(CP*CE*(TE*CD-SA*SD)+CA*SD*SP)

        devRA.append(dalpha*R2A)
        devDEC.append(ddec*R2A)

        inp[ii]['RightAscension'] = inp[ii]['RightAscension'] + dalpha
        inp[ii]['Declination'] = inp[ii]['Declination'] + ddec

        # TEST CASE
        # alpha = DTR*41.5472
        # dec = DTR*49.3485
        # T = 0.2886705
        # print('eps =', R2D*eps - 23.436)
        # print('lonS =',R2D*lonS - 231.328)
        # print('ecc =',e - 0.01669649)
        # print('pi =',R2D*pi - 103.434)
        # print('alpha =',R2A*dalpha - 30.045)
        # print('dec =',R2A*ddec - 6.697)

    with open(obs_data, 'w') as fp:
        json.dump(inp, fp, indent=4)

    if plot:
        TITLE = 'Stellar Aberration Correction'
        # Plot Filter Measurement Residuals
        title1 = 'Right Ascension'
        trace1 = go.Scatter(x=tstr, y=devRA, mode='markers', marker=dict(
            color='rgb(51, 153, 255)'), showlegend=False)
        title2 = 'Declination'
        trace2 = go.Scatter(x=tstr, y=devDEC, mode='markers', marker=dict(
            color='rgb(0, 204, 0)'), showlegend=False)
        # Create subplot
        fig = tools.make_subplots(rows=2, cols=1, subplot_titles=(
            title1, title2), shared_xaxes=True)

        # Assign traces to subplots
        fig.append_trace(trace1, 1, 1)
        fig.append_trace(trace2, 2, 1)
        # Configure Layout
        fig['layout'].update(title=TITLE)
        fig['layout']['xaxis1'].update(title='Time')
        fig['layout']['yaxis1'].update(title='Aberration Effect [arcsec]')
        fig['layout']['yaxis2'].update(title='Aberration Effect [arcsec]')
        # Execute
        plotly.offline.plot(fig, filename='Plots/abereffect.html', auto_open=True)

    return


def plotEstSensorError(od_cfg, obs_data):
    """Plot residuals between observation file and reference data in measurement space
    :param od_cfg: od_cfg filename
    :type obs_data: obs_data filename
    """

    with open(od_cfg, "r") as fp:
            cfg = json.load(fp)
    with open(obs_data, "r") as fp:
            inp = json.load(fp)

    key = tuple(cfg["Measurements"].keys())

    # Get site location
    key2 = tuple(cfg["Stations"].keys())
    lat = cfg['Stations'][key2[0]]['Latitude']
    lon = cfg['Stations'][key2[0]]['Longitude']
    alt = cfg['Stations'][key2[0]]['Altitude']

    tstr, pvRef, ang1, ang2 = [], [], [], []
    for i in inp:
            # Get time string
            tstr.append(i["Time"])
            # Get reference state
            PV = np.array(i["TrueState"]["Cartesian"][:6])
            pvRef.append(PV)
            # Get obs
            ang1.append(i[key[0]])
            ang2.append(i[key[1]])

    obs = np.array([ang1, ang2])
    pvRef = np.array(pvRef)
    numObs = np.size(ang1)

    smoothMeas, smoothRes, refMeas, refRes = [], [], [], []
    sig1 = np.array(cfg['Measurements'][key[0]]['Error'])
    sig2 = np.array(cfg['Measurements'][key[1]]['Error'])
    for ii in range(0, numObs, 1):
        ti = tstr[ii]
        sigma = [sig1, sig2]

        pv = pvRef[ii, :].tolist()
        angular = obs[:, ii].tolist()
        if key[0] == 'RightAscension':
            output = pv2radec(od_cfg, lat, lon, alt, ti, angular, sigma, pv)
        if key[0] == 'Azimuth':
            output = pv2azel(od_cfg, lat, lon, alt, ti, angular, sigma, pv)
        if key[0] == 'Range':
            range1 = obs[0, ii].tolist()
            rangerate1 = obs[1, ii].tolist()
            output1 = pv2range(od_cfg, lat, lon, alt, ti, range1, sig1, pv)
            output2 = pv2rangerate(od_cfg, lat, lon, alt,
                                   ti, rangerate1, sig2, pv)
            output = [output1, output2]
        refMeas.append(output)
        output = np.squeeze(np.array(output))
        angular = np.array(angular)
        refRes.append(angular - output)

    refRes = np.array(refRes)
    meanstdANGref = meanstdANG(refRes)

    angles = ("Azimuth", "Elevation", "RightAscension", "Declination")
    if (key[0] in angles and key[1] in angles):
        refRes = refRes*648000.0/math.pi
        meanstdANGref = meanstdANG(refRes)
        units = ("[arcsec]", "[arcsec]")
    else:
        if ("PositionVelocity" in key):
            
            units = ("[m]", "[m]", "[m]", "[m/s]", "[m/s]", "[m/s]")
        else:
            units = ("[m]", "[m/s]")

    TITLE = 'Sensor Measurement - Reference Measurement Residuals'
    # Plot Filter Measurement Residuals
    title1 = key[0] + ' Error: Mean = ' + \
        format(meanstdANGref[0], '.4f') + ' | STD = ' + \
        format(meanstdANGref[1], '.4f')
    trace1 = go.Scatter(x=tstr, y=refRes[:, 0], mode='markers', name=key[0], marker=dict(
        color='rgb(51, 153, 255)'), showlegend=False)
    title2 = key[1] + ' Declination Error: Mean = ' + \
        format(meanstdANGref[2], '.4f') + ' | STD = ' + \
        format(meanstdANGref[3], '.4f')
    trace2 = go.Scatter(x=tstr, y=refRes[:, 1], mode='markers', name=key[1], marker=dict(
        color='rgb(0, 204, 0)'), showlegend=False)
    # Create subplot
    fig = tools.make_subplots(rows=2, cols=1, subplot_titles=(
        title1, title2), shared_xaxes=True)
    # Assign traces to subplots
    fig.append_trace(trace1, 1, 1)
    fig.append_trace(trace2, 2, 1)
    # Configure Layout
    fig['layout'].update(title=TITLE)
    fig['layout']['xaxis1'].update(title='Time')
    fig['layout']['yaxis1'].update(title=key[0]+' '+units[0])
    fig['layout']['yaxis2'].update(title=key[1]+' '+units[1])
    # Execute
    plotly.offline.plot(fig, filename='Plots/measerror.html', auto_open=True)

    return


def parseAstroAndRefData(fnAstro, fnRef, fnObs, fnODout, tBias):

    fnODin = 'data/radec_od_cfg.json'  # Use as generic cfg to edit
    # Load astro formatted .mat file into python
    mat_contents = sio.loadmat(fnAstro, struct_as_record=False)
    astro = mat_contents['astro']
    astro = astro[0, 0]
    # Load measurements, time, and site location
    RA = astro.ra_stellar
    DEC = astro.dec_stellar
    numObs = np.size(RA)
    obs_mjd = astro.mjd
    obs_Time = Time(obs_mjd, format='mjd', scale='utc') - timedelta(seconds = tBias)
    obs_isot = np.squeeze(obs_Time.isot)
    obs_mjd = obs_Time.mjd
    # Get lat, lon, alt, for station location in od_cfg (radians and meters)
    lat = float(astro.site_lat*math.pi/180)
    lon = float(astro.site_lon*math.pi/180)
    alt = float(astro.site_alt*1000)
    # Load reference satellite data for initial position, velocity, and time for od_cfg (J200 meters)
    filestore = sio.loadmat(fnRef)
    refSat = filestore['x_refsat']
    PV = refSat[:, 2:8]*1000
    ref_mjd = refSat[:, 0] - 2400000.5
    PVobs = interpolationPV(np.squeeze(obs_mjd), ref_mjd, PV, 10)
    # DATA SEGMENTATION
    startMeasurement = int(0)
    endMeasurement = int(numObs)
    startTime = obs_isot[startMeasurement]
    endTime = obs_isot[endMeasurement-1]
    writeObsJSON(fnObs, startMeasurement, endMeasurement,
                 obs_isot, PVobs, RA, DEC)
    writeODcfgJSON(fnODin, fnODout, startTime, endTime,
                   startMeasurement, PVobs, lat, lon, alt)


def lla2eci(gslat, gslon, gsalt, time):
    """Converts a lat, lon, alt to eciJ200 coordinates
    :param gslat: geodetic latitude of station [rad]
    :param gslon: geodetic longitude of station [rad]
    :param gsalt: altitude of station [m]
    :param time: time stamp of observation
    :type gslat: float
    :type gslon: float
    :type gsalt: float
    :type time: str
    :return eciposvel: ECI J2000 position and velocity of station
    :rtype: list
    """
    ecij2000 = lla2eciJ2000(gslat, gslon, gsalt, time)
    return(ecij2000)


def estTimeBiasRADEC(obs_data,od_cfg):
    """ 
    estimates a time bias using reference data for right ascension and declination measurements
    """
    with open(od_cfg, "r") as fp:
        cfg = json.load(fp)
    with open(obs_data, "r") as fp:
        inp = json.load(fp)
    key = tuple(cfg["Measurements"].keys())
    # Get site location
    key2 = tuple(cfg["Stations"].keys())
    lat = cfg['Stations'][key2[0]]['Latitude']
    lon = cfg['Stations'][key2[0]]['Longitude']
    alt = cfg['Stations'][key2[0]]['Altitude']
    sig1 = np.array(cfg['Measurements'][key[0]]['Error'])
    sig2 = np.array(cfg['Measurements'][key[1]]['Error'])
    tstr, pvRef, ang1, ang2 = [], [], [], []
    for i in inp:
            # Get time string
            tstr.append(i["Time"])
            # Get reference state
            PV = np.array(i["TrueState"]["Cartesian"][:6])
            pvRef.append(PV)
            # Get obs
            ang1.append(i[key[0]])
            ang2.append(i[key[1]])
    # Raw observation at time (ti + tb)
    obs = np.array([ang1, ang2])
    # Reference sat position and velocity at time (ti + tb)
    pvRef = np.array(pvRef)
    # Number of observations
    numObs = np.size(ang1)
    a0 = 0
    a1 = 0
    a2 = 0
    rms = 0
    rms_prev = 999
    tb_eps = 1e-10
    count = 0
    maxcount = 10
    timeObj = Time(tstr, format='isot', scale='utc')
    obs_mjd = timeObj.mjd
    tObs = (obs_mjd - obs_mjd[0]) * 86400
    degree = 1
    while abs(rms-rms_prev) > tb_eps and count < maxcount : 
        z = np.zeros([2*numObs,1])
        H = np.zeros([2*numObs,degree])
        W = np.identity(2*numObs) * 1/(math.pi * 5 / 3600 / 180)**2
        count = count + 1
        for ii in range(0, numObs, 1):
            tBias = a0 + a1*tObs[ii] + a2*tObs[ii]**2
            t = tObs[ii] - tBias
            sigma = [sig1, sig2]
            pv = interpolationPV(t, tObs, pvRef, 10)
            pv = np.squeeze(pv)
            pv = pvRef[ii,:]
            obs_Time = Time(obs_mjd[ii], format='mjd', scale='utc') - timedelta(seconds=tBias)
            obs_isot = obs_Time.isot
            refObs = pv2radec(od_cfg, lat, lon, alt, obs_isot, obs[:, ii].tolist(), sigma, pv.tolist())
            epsRA = obs[0,ii] - refObs[0]
            epsDEC = obs[1,ii] - refObs[1]
            z[2*ii-1] = epsRA
            z[2*ii] = epsDEC
            # print(epsRA*math.pi*648000.0/math.pi)
            vx = pv[3]
            vy = pv[4]
            vz = pv[5]
            vVec = pv[3:]
            rSite = lla2eci(lat, lon, alt, obs_isot)
            rhoVec = pv[:3] - rSite - pv[3:]*tBias
            rhox = rhoVec[0]
            rhoy = rhoVec[1]
            rhoz = rhoVec[2]
            rho = np.linalg.norm(rhoVec[:3])
            rhoHat = rhoVec/rho
            for jj in range(0,degree,1):
                depsRAdaj = -( rhox*vy - rhoy*vx ) / (rhox**2 + rhoy**2) * t**jj
                H[2*ii-1,jj] = depsRAdaj
                depsDECdaj = -( rho*vz - rhoz*np.dot(rhoHat,vVec) ) / (rho**2*math.sqrt(1-rhoHat[2]**2)) * t**jj
                H[2*ii,jj] = depsDECdaj
        HtH = H.T @ W @ H
        Phat = np.linalg.inv(HtH)
        xHat = Phat @ H.T @ W @ z
        a0 = float(xHat[0])
        if degree > 1:
            a1 = float(xHat[1])
        if degree > 2:
            a2 = float(xHat[2])
        rms_prev = rms
        rms = np.sqrt( sum( z**2, 1 ) / 2 / numObs )
    return a0

def addTimeBias(obs_data, od_cfg, tBias):
    """ applies a given time bias to the obs data and configure file """
    with open(od_cfg, "r") as fp:
        cfg = json.load(fp)
    with open(obs_data, "r") as fp:
        inp = json.load(fp)
    pvRef, refMJD, obsMJD = [], [], []
    for i in inp:
        refTimeString = i['Time']
        refTime = Time(refTimeString, format='isot', scale='utc')
        refMJD.append(refTime.mjd)
        obsTime = Time(refTimeString, format='isot', scale='utc') - timedelta(seconds = tBias)
        obsMJD.append(obsTime.mjd)
        i['Time'] = obsTime.isot + 'Z'
        pvRef.append(i['TrueState']['Cartesian'])
    pvRef = np.array(pvRef)
    refMJD = np.array(refMJD)
    startTime = inp[0]['Time']
    endTime = inp[-1]['Time']
    cfg['Propagation']['Start'] = startTime
    cfg['Propagation']['End'] = endTime
    PVobs = interpolationPV(np.squeeze(obsMJD), refMJD, pvRef, 10)
    count = 0
    for i in inp:
        i['TrueState']['Cartesian'] = PVobs[count,:].tolist()
        count = count + 1
    with open(obs_data, 'w') as fp:
        json.dump(inp, fp, indent=4)
    with open(od_cfg, 'w') as fp:
        json.dump(cfg, fp, indent=4)
    return


def estTimeBiasRange(obs_data,od_cfg):
    """ NOT COMPLETE, estimate time bias using reference data for range, az, and el measurements"""
    with open(od_cfg, "r") as fp:
        cfg = json.load(fp)
    with open(obs_data, "r") as fp:
        inp = json.load(fp)
    key = tuple(cfg["Measurements"].keys())
    # Get site location
    key2 = tuple(cfg["Stations"].keys())
    lat = cfg['Stations'][key2[0]]['Latitude']
    lon = cfg['Stations'][key2[0]]['Longitude']
    alt = cfg['Stations'][key2[0]]['Altitude']
    sig1 = np.array(cfg['Measurements'][key[0]]['Error'])
    sig2 = np.array(cfg['Measurements'][key[1]]['Error'])
    tstr, pvRef, ang1, ang2 = [], [], [], []
    for i in inp:
            # Get time string
            tstr.append(i["Time"])
            # Get reference state
            PV = np.array(i["TrueState"]["Cartesian"][:6])
            pvRef.append(PV)
            # Get obs
            ang1.append(i[key[0]])
            ang2.append(i[key[1]])
    # Raw observation at time (ti + tb)
    obs = np.array([ang1, ang2])
    # Reference sat position and velocity at time (ti + tb)
    pvRef = np.array(pvRef)
    # Number of observations
    numObs = np.size(ang1)
    a0 = 0
    a1 = 0
    a2 = 0
    rms = 0
    rms_prev = 999
    tb_eps = 1e-10
    count = 0
    maxcount = 10
    timeObj = Time(tstr, format='isot', scale='utc')
    obs_mjd = timeObj.mjd
    tObs = (obs_mjd - obs_mjd[0]) * 86400
    degree = 1
    while abs(rms-rms_prev) > tb_eps and count < maxcount : 
        z = np.zeros([2*numObs,1])
        H = np.zeros([2*numObs,degree])
        W = np.identity(2*numObs) * 1/(math.pi * 5 / 3600 / 180)**2
        count = count + 1
        for ii in range(0, numObs, 1):
            tBias = a0 + a1*tObs[ii] + a2*tObs[ii]**2
            t = tObs[ii] - tBias
            sigma = [sig1, sig2]
            pv = interpolationPV(t, tObs, pvRef, 10)
            pv = np.squeeze(pv)
            pv = pvRef[ii,:]
            obs_Time = Time(obs_mjd[ii], format='mjd', scale='utc') - timedelta(seconds=tBias)
            obs_isot = obs_Time.isot
            refObs = pv2radec(od_cfg, lat, lon, alt, obs_isot, obs[:, ii].tolist(), sigma, pv.tolist())
            refObs1 = pv2range(od_cfg, lat, lon, alt, obs_isot, obs[0, ii].tolist(), sig1, pv.tolist())
            refObs2 = pv2rangerate(od_cfg, lat, lon, alt, obs_isot, obs[1, ii].tolist(), sig2, pv.tolist())
            epsRA = obs[0,ii] - refObs1
            epsDEC = obs[1,ii] - refObs
            z[2*ii-1] = epsRA
            z[2*ii] = epsDEC
            # print(epsRA*math.pi*648000.0/math.pi)
            vx = pv[3]
            vy = pv[4]
            vz = pv[5]
            vVec = pv[3:]
            rSite = lla2eci(lat, lon, alt, obs_isot)
            rhoVec = pv[:3] - rSite - pv[3:]*tBias
            rhox = rhoVec[0]
            rhoy = rhoVec[1]
            rhoz = rhoVec[2]
            rho = np.linalg.norm(rhoVec[:3])
            rhoHat = rhoVec/rho
            for jj in range(0,degree,1):
                depsRAdaj = -( rhox*vy - rhoy*vx ) / (rhox**2 + rhoy**2) * t**jj
                H[2*ii-1,jj] = depsRAdaj
                depsDECdaj = -( rho*vz - rhoz*np.dot(rhoHat,vVec) ) / (rho**2*math.sqrt(1-rhoHat[2]**2)) * t**jj
                H[2*ii,jj] = depsDECdaj
        HtH = H.T @ W @ H
        Phat = np.linalg.inv(HtH)
        xHat = Phat @ H.T @ W @ z
        a0 = float(xHat[0])
        if degree > 1:
            a1 = float(xHat[1])
        if degree > 2:
            a2 = float(xHat[2])
        rms_prev = rms
        rms = np.sqrt( sum( z**2, 1 ) / 2 / numObs )

    return a0