###
# Plot spin with respect to time
###

import numpy as np
import matplotlib.pyplot as plt

yr = 365.25

spinfile = open("spin.out","r")
data = np.loadtxt(spinfile)
spinfile.close()

time = data[:,0]
spindata= data[:,1:]

nt = len(time)
nbod = int(spindata.shape[1]/3)

spins = np.zeros((nt,nbod))
for t in range(nt):
    for j in range(nbod):
        k = j*3
        sx = spindata[t, k]
        sy = spindata[t, k+1]
        sz = spindata[t, k+2]
        spins[t, j] = np.sqrt(sx*sx + sy*sy + sz*sz)

nbodspin = np.sum(spins[0]>0)

print(f"{nbod} massive bodies, {nbodspin} with a spin")

for k in range(nbod):
    if spins[0,k] > 0:
        rotationperiod = (2*np.pi/spins[:,k])*yr
        if nbodspin > 1:
            plt.subplot(1,nbodspin,k+1)
        plt.plot(time, rotationperiod, '.')
        plt.xlabel("Time (yr)")
        plt.ylabel("Rotation period (d)")
        plt.title(f"Body #{k+1}")

plt.tight_layout()
plt.show()
