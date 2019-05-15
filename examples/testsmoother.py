import sys
import time
import smoother
import matplotlib.pyplot as plt

print("Smoother start : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))

filename = sys.argv[1]

xPre, xPost, PPre, PPost, filetype = smoother.dataConfigJSON(filename)
xS, PS, T = smoother.rtsSmoother(xPre, xPost, PPre, PPost, filetype)

residuals = abs(xPost[0:3, :] - xS[0:3, :])
plt.plot(T, residuals[0, :], label='x')
plt.plot(T, residuals[1, :], label='y')
plt.plot(T, residuals[2, :], label='z')
plt.xlabel('Measurement Number')
plt.ylabel('Residuals (postFit - smoothed)')
plt.title('Smoothing Results')
plt.legend()

print("Smoother end   : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))

plt.show()


