from matplotlib import pyplot as plt
import numpy as np


def Vpot(x):
    return (x**2 - 2.5) * x**2
    # return 0.5*x**2


a = 10
N = 1000
mu = 2
sigma = 1
x = np.linspace(-a / 2, a / 2, N)
psi = np.exp(-(x - mu)**2 / (2 * sigma)) + np.exp(-(x + mu)**2 / (2 * sigma))
psi *= 500
p = psi * psi
V = Vpot(x)
plt.plot(x, psi, label='p')
plt.plot(x, V, label='v')
plt.legend()
plt.grid()
plt.show()
