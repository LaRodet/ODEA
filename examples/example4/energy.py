###
# Plot energy with respect to time
###

import numpy as np
import matplotlib.pyplot as plt

energyfile = open("energy.out","r")
data = np.loadtxt(energyfile)
energyfile.close()

time = data[:,0]
energy = data[:,1]

plt.plot(time,energy/energy[0],'.')
plt.xlabel("Time (yr)")
plt.ylabel("Energy (normalized)")

plt.show()
