from matplotlib import pyplot as plt
import numpy as np

datafile = "./data/stats_1.2.1.dat"
sums1, sums2, sums10, sums100 = np.loadtxt(
    datafile, usecols=(0, 1, 2, 3), delimiter=' ', unpack='true')

fig, axes = plt.subplots(2, 2, figsize=(7, 5))
plt.tight_layout()
fig.suptitle("Sums of Uniformly Distributed Random Numbers")

axes[0, 0].hist(sums1, bins=100, ec='black', density=True)
axes[0, 0].title.set_text("N = 1")

axes[0, 1].hist(sums2, bins=100, ec='black', density=True)
axes[0, 1].title.set_text("N = 2")

axes[1, 0].hist(sums10, bins=100, ec='black', density=True)
axes[1, 0].title.set_text("N = 10")

axes[1, 1].hist(sums100, bins=100, ec='black', density=True)
axes[1, 1].title.set_text("N = 100")

plt.tight_layout()
plt.show()
