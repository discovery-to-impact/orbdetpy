import numpy as np
import json

def dataConfigJSON(filename):
    filetype = 'json'
    with open(filename) as f:
        jsonData = json.load(f)

    alpha = jsonData["Measurement"]["alpha"]
    dec = jsonData["Measurement"]["dec"]
    alphaDot = jsonData["Measurement"]["alphaDot"]
    decDot = jsonData["Measurement"]["decDot"]
    zVec = np.array([[alpha, dec, alphaDot, decDot]]).T

    rx = jsonData["Station"]["rx"]
    ry = jsonData["Station"]["ry"]
    rz = jsonData["Station"]["rz"]
    vx = jsonData["Station"]["vx"]
    vy = jsonData["Station"]["vy"]
    vz = jsonData["Station"]["vz"]
    rVecS = np.array([[rx, ry, rz]]).T
    rDotVecS = np.array([[vx, vy, vz]]).T
    aMin = jsonData["Constraints"]["aMin"]
    aMax = jsonData["Constraints"]["aMax"]
    eMax = jsonData["Constraints"]["eMax"]

    sigmaTildeRho = jsonData["Constraints"]["sigmaTildeRho"]
    sigmaTildeRhoDot = jsonData["Constraints"]["sigmaTildeRhoDot"]
    separationDistance = jsonData["Constraints"]["separationDistance"]

    return zVec, rVecS, rDotVecS, aMin, aMax, eMax, sigmaTildeRho, sigmaTildeRhoDot, separationDistance

def constrainSMA(rhoSet, zVec, rVecS, rDotVecS, energy):
    alpha = zVec[0]
    dec = zVec[1]
    alphaDot = zVec[2]
    decDot = zVec[3]

    first = 1

    for ii in range(0, np.size(rhoSet), 1):
        rho = rhoSet[ii]

        uVecRho = [np.cos(alpha) * np.cos(dec), np.sin(alpha) * np.cos(dec), np.sin(dec)]
        uVecAlpha = [-np.sin(alpha) * np.cos(dec), np.cos(alpha) * np.cos(dec), np.array([0])]
        uVecDec = [-np.cos(alpha) * np.sin(dec), -np.sin(alpha) * np.sin(dec), np.cos(dec)]

        uVecRho = np.array(uVecRho)
        uVecAlpha = np.array(uVecAlpha)
        uVecDec = np.array(uVecDec)

        qVec = rVecS.T
        qDotVec = rDotVecS.T

        w0 = np.linalg.norm(qVec) ** 2
        w1 = 2 * np.vdot(qDotVec, uVecRho)
        w2 = alphaDot ** 2 * np.cos(dec) ** 2 + decDot ** 2
        w3 = 2 * alphaDot * np.vdot(qDotVec, uVecAlpha) + 2 * decDot * np.vdot(qDotVec, uVecDec)
        w4 = np.linalg.norm(qDotVec) ** 2
        w5 = 2 * np.vdot(qVec, uVecRho)

        mu = 398600.4415  # [km^3/s^2]

        Frho = w2 * rho ** 2 + w3 * rho + w4 - 2 * mu / np.sqrt(rho ** 2 + w5 * rho + w0)
        rhoDotPos = -w1 / 2 + np.sqrt((w1 / 2) ** 2 - Frho + 2 * energy)
        rhoDotNeg = -w1 / 2 - np.sqrt((w1 / 2) ** 2 - Frho + 2 * energy)

        if not np.isnan(rhoDotPos) and not np.isnan(rhoDotNeg):
            if first == 1:
                rhoDot = np.array([[rhoDotPos, rhoDotNeg]])
                rhoStore = np.array([[rhoSet[ii]]])
                first = 0
            else:
                rhoDot = np.concatenate((rhoDot, np.array([[rhoDotPos, rhoDotNeg]])), axis=0)
                rhoStore = np.concatenate((rhoStore, np.array([[rhoSet[ii]]])), axis=0)

    return rhoStore, np.squeeze(rhoDot, axis=2)


def constrainECC(rhoSet, zVec, rVecS, rDotVecS, e):
    alpha = zVec[0]
    dec = zVec[1]
    alphaDot = zVec[2]
    decDot = zVec[3]

    first = 1

    for ii in range(0, np.size(rhoSet), 1):
        rho = rhoSet[ii]
        uVecRho = [np.cos(alpha) * np.cos(dec), np.sin(alpha) * np.cos(dec), np.sin(dec)]
        uVecAlpha = [-np.sin(alpha) * np.cos(dec), np.cos(alpha) * np.cos(dec), np.array([0])]
        uVecDec = [-np.cos(alpha) * np.sin(dec), -np.sin(alpha) * np.sin(dec), np.cos(dec)]

        uVecRho = np.array(uVecRho)
        uVecAlpha = np.array(uVecAlpha)
        uVecDec = np.array(uVecDec)

        qVec = rVecS
        qDotVec = rDotVecS

        w0 = np.linalg.norm(qVec) ** 2
        w1 = 2 * np.vdot(qDotVec, uVecRho)
        w2 = alphaDot ** 2 * np.cos(dec) ** 2 + decDot ** 2
        w3 = 2 * alphaDot * np.vdot(qDotVec, uVecAlpha) + 2 * decDot * np.vdot(qDotVec, uVecDec)
        w4 = np.linalg.norm(qDotVec) ** 2
        w5 = 2 * np.vdot(qVec, uVecRho)

        mu = 398600.4415  # [km^3/s^2]

        Frho = w2 * rho ** 2 + w3 * rho + w4 - 2 * mu / np.sqrt(rho ** 2 + w5 * rho + w0)

        hVec1 = np.cross(qVec, uVecRho, axis=0)
        hVec2 = np.cross(uVecRho, (alphaDot * uVecAlpha + decDot * uVecDec), axis=0)
        hVec3 = np.cross(uVecRho, qDotVec, axis=0) + np.cross(qVec, (alphaDot * uVecAlpha + decDot * uVecDec), axis=0)
        hVec4 = np.cross(qVec, qDotVec, axis=0)

        c0 = np.linalg.norm(hVec1) ** 2
        c1 = 2 * np.vdot(hVec1, hVec2)
        c2 = 2 * np.vdot(hVec1, hVec3)
        c3 = 2 * np.vdot(hVec1, hVec4)
        c4 = np.linalg.norm(hVec2) ** 2
        c5 = 2 * np.vdot(hVec2, hVec3)
        c6 = 2 * np.vdot(hVec2, hVec4) + np.linalg.norm(hVec3) ** 2
        c7 = 2 * np.vdot(hVec3, hVec4)
        c8 = np.linalg.norm(hVec4) ** 2

        Prho = c1 * rho ** 2 + c2 * rho + c3
        Urho = c4 * rho ** 4 + c5 * rho ** 3 + c6 * rho ** 2 + c7 * rho + c8

        a0 = Frho * Urho + mu ** 2 * (1 - e ** 2)
        a1 = Frho * Prho + w1 * Urho
        a2 = Urho + c0 * Frho + w1 * Prho
        a3 = Prho + c0 * w1
        a4 = c0

        eqnCoff = [a4, a3, a2, a1, a0]
        rhoDotSolve = np.roots(eqnCoff)

        if np.isreal(rhoDotSolve[2]) and np.isreal(rhoDotSolve[3]):
            if first == 1:
                rhoDot = np.array([[rhoDotSolve[2], rhoDotSolve[3]]])
                rhoStore = np.array([[rhoSet[ii]]])
                first = 0
            else:
                rhoDot = np.concatenate((rhoDot, np.array([[rhoDotSolve[2], rhoDotSolve[3]]])), axis=0)
                rhoStore = np.concatenate((rhoStore, np.array([[rhoSet[ii]]])), axis=0)

    rhoDotOutput = np.real(rhoDot)
    rhoDotPos = rhoDotOutput[rhoDotOutput >= 0]
    rhoDotNeg = rhoDotOutput[rhoDotOutput < 0]
    rhoDotOut = np.array([rhoDotPos.T, rhoDotNeg.T])
    return rhoStore, rhoDotOut.T


def truncateCAR(rhoSet, rhoSet1, rhoSet2, rhoSet3, rhoDotSMA_Max, rhoDotSMA_Min, rhoDotECC):

    R = 6378.1363

    rangeRateUpper = np.array([])
    rangeRateLower = np.array([])

    leftRho = np.amax([np.amax(rhoSet2), np.amin(rhoSet3)])
    rightRho = np.amin([np.amax(rhoSet1), np.amax(rhoSet3)])
    leftIndex = np.where(leftRho == rhoSet)
    rightIndex = np.where(rightRho == rhoSet)
    leftIndex = np.ndarray.item(leftIndex[0])
    rightIndex = np.ndarray.item(rightIndex[0])
    rhoSetTrunc = rhoSet[leftIndex:rightIndex+1]

    # Get max range rate value given range
    for ii in range(0, np.size(rhoSetTrunc), 1):

        set1RangeRate = np.array([])
        set2RangeRate = np.array([])
        set3RangeRate = np.array([])

        index1 = np.where(rhoSetTrunc[ii] == rhoSet1)
        index2 = np.where(rhoSetTrunc[ii] == rhoSet2)
        index3 = np.where(rhoSetTrunc[ii] == rhoSet3)

        if index1[0].size > 0:
            set1RangeRate = rhoDotSMA_Max[index1[0], 0]

        if index2[0].size > 0:
            set2RangeRate = rhoDotSMA_Min[index2[0], 0]

        if index3[0].size > 0:
            set3RangeRate = np.array([np.amax(rhoDotECC[index3[0], :])])

        rangeRateMin = min([set1RangeRate, set2RangeRate, set3RangeRate])

        if rangeRateUpper.size == 0:
            rangeRateUpper = rangeRateMin
        else:
            if index3[0] == 520:
                a = 0
            rangeRateUpper = np.concatenate((rangeRateUpper, rangeRateMin), axis=0)

    for ii in range(0, np.size(rhoSetTrunc), 1):

        set1RangeRate = np.array([])
        set2RangeRate = np.array([])
        set3RangeRate = np.array([])

        index1 = np.where(rhoSetTrunc[ii] == rhoSet1)
        index2 = np.where(rhoSetTrunc[ii] == rhoSet2)
        index3 = np.where(rhoSetTrunc[ii] == rhoSet3)

        if index1[0].size > 0:
            set1RangeRate = rhoDotSMA_Max[index1[0], 1]

        if index2[0].size > 0:
            set2RangeRate = rhoDotSMA_Min[index2[0], 1]

        if index3[0].size > 0:
            set3RangeRate = np.array([np.amin(rhoDotECC[index3[0], :])])

        if rangeRateLower.size == 0:
            rangeRateMin = np.amax([set1RangeRate,set2RangeRate,set3RangeRate])
            rangeRateLower = rangeRateMin
        else:
            rangeRateMin = max([set1RangeRate, set2RangeRate, set3RangeRate])
            rangeRateLower = np.concatenate((rangeRateLower, rangeRateMin), axis=0)

    rhoDotTrunc = np.array([rangeRateUpper.T,rangeRateLower.T])
    rhoSetTrunc = np.array([rhoSet[leftIndex:rightIndex+1]])

    return rhoSetTrunc.T, rhoDotTrunc.T
