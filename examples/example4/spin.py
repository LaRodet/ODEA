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
sx = np.zeros((nt,nbod))
sy = np.zeros((nt,nbod))
sz = np.zeros((nt,nbod))
for t in range(nt):
    for j in range(nbod):
        sx[t,j] = spindata[t, j]
        sy[t,j] = spindata[t, 3+j]
        sz[t,j] = spindata[t, 6+j]
        spins[t, j] = np.sqrt(sx[t,j]**2 + sy[t,j]**2 + sz[t,j]**2)

### Rotation rate

nbodspin = np.sum(spins[0]>0)

print(f"{nbod} massive bodies, {nbodspin} with a spin")

for k in range(nbod):
    if spins[0,k] > 0:
        rotationperiod = (2*np.pi/spins[:,k])*yr
        if nbodspin > 1:
            plt.subplot(1,nbodspin,k+1)
        plt.plot(time, rotationperiod, '.-')
        plt.xlabel("Time (yr)")
        plt.ylabel("Rotation period (d)")
        plt.title(f"Body #{k+1}")

plt.tight_layout()

### sx

plt.figure()

nbodspin = np.sum(np.array([np.any(sx[:,j]>0) for j in range(nbod)]))

for k in range(nbod):
    if np.any(sx[:,k]) > 0:
        if nbodspin > 1:
            plt.subplot(1,nbodspin,k+1)
        plt.plot(time, sx[:,k]/yr, '.-')
        plt.xlabel("Time (yr)")
        plt.ylabel(r"Spin x coordinate (d$^{-1}$)")
        plt.title(f"Body #{k+1}")

plt.tight_layout()

### sy

plt.figure()

nbodspin = np.sum(np.array([np.any(sy[:,j]>0) for j in range(nbod)]))

for k in range(nbod):
    if np.any(sy[:,k]) > 0:
        if nbodspin > 1:
            plt.subplot(1,nbodspin,k+1)
        plt.plot(time, sy[:,k]/yr, '.-')
        plt.xlabel("Time (yr)")
        plt.ylabel(r"Spin y coordinate (d$^{-1}$)")
        plt.title(f"Body #{k+1}")

plt.tight_layout()

### sz

plt.figure()

nbodspin = np.sum(np.array([np.any(sz[:,j]>0) for j in range(nbod)]))

for k in range(nbod):
    if np.any(sz[:,k]) > 0:
        if nbodspin > 1:
            plt.subplot(1,nbodspin,k+1)
        plt.plot(time, sz[:,k]/yr, '.-')
        plt.xlabel("Time (yr)")
        plt.ylabel(r"Spin z coordinate (d$^{-1}$)")
        plt.title(f"Body #{k+1}")

plt.tight_layout()


plt.show()
