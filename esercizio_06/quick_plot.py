from matplotlib import pyplot as plt
import numpy as np

datafile = "./output.heat.0"
val = np.loadtxt(datafile, usecols=(1), unpack='true')
x = np.arange(val.size)
plt.figure(figsize=(10, 5))
plt.grid(True)
plt.plot(x, val)
plt.show()

avg, err = np.loadtxt(datafile, usecols=(2, 3), unpack='true')

x = np.arange(avg.size)
plt.figure(figsize=(10, 5))
plt.grid(True)
plt.errorbar(x, avg, yerr=err)
plt.show()
