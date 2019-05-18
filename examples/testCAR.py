import sys
import time
import admissibleRegion as AR
import matplotlib.pyplot as plt

print("CAR start : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))

filename = sys.argv[1]

[zVec, rVecS, rDotVecS, aMin, aMax, eMax, sigmaTildeRho, sigmaTildeRhoDot, separationDistance] = AR.dataConfigJSON(filename)

# Constants used
R = 6378.1363  # [km]
mu = 398600.4415  # [km^3/s^2]

rangeMin = aMin - R
rangeMax = aMax + R
rhoSet = range(0, int(7 * R), separationDistance)
energyMax = -mu / (2 * aMax)
energyMin = -mu / (2 * aMin)

# SMA Constraint Max
[rhoSet1, rhoDotSMA_Max] = AR.constrainSMA(rhoSet, zVec, rVecS, rDotVecS, energyMax)
# SMA Constraint Min
[rhoSet2, rhoDotSMA_Min] = AR.constrainSMA(rhoSet, zVec, rVecS, rDotVecS, energyMin)
# Eccentricity Constraint
[rhoSet3, rhoDotECC] = AR.constrainECC(rhoSet, zVec, rVecS, rDotVecS, eMax)
# Truncate
[rhoSetTrunc,rhoDotTrunc] = AR.truncateCAR(rhoSet, rhoSet1, rhoSet2, rhoSet3, rhoDotSMA_Max, rhoDotSMA_Min, rhoDotECC)

plt.figure()
plt.plot(rhoSet1 / R, rhoDotSMA_Max[:, 0], 'blue', label="$a_{max} = %i km$" % aMax)
plt.plot(rhoSet1 / R, rhoDotSMA_Max[:, 1], 'blue')
plt.plot(rhoSet2 / R, rhoDotSMA_Min[:, 0], 'orange', label="$a_{min} = %i km$" % aMin)
plt.plot(rhoSet2 / R, rhoDotSMA_Min[:, 1], 'orange')
plt.plot(rhoSet3 / R, rhoDotECC[:, 0], 'green', label="$e_{max} = %.1f$" % eMax)
plt.plot(rhoSet3 / R, rhoDotECC[:, 1], 'green')
plt.xlabel('Range [Earth Radius]')
plt.ylabel("Range Rate [$\dfrac{km}{s}$]")
plt.title('Semi-Major Axis and Eccentricity Constrained Admissible Region')
plt.grid()
plt.legend()

plt.figure()
plt.plot(rhoSetTrunc / R, rhoDotTrunc[:,0], 'green')
plt.plot(rhoSetTrunc / R, rhoDotTrunc[:,1], 'green')
plt.xlabel('Range [Earth Radius]')
plt.ylabel("Range Rate [$\dfrac{km}{s}$]")
plt.title('Semi-Major Axis and Eccentricity Constrained Admissible Region Truncated')
plt.grid()

print("CAR end   : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))

plt.show()