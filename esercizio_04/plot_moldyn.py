import numpy
import seaborn
import matplotlib.pyplot

sigma = 3.4 * 10 ** (-10)  # m
k_B = 1.380649 * 10 ** (-23)  # J/K
T = 120
epsilon = T * k_B  # J
m = 39.948 * 1.66 * 10 ** (-27)  # mass of argon, kg
t = sigma * numpy.sqrt(m / epsilon)
p = epsilon / sigma ** 3  # Pa

lj_scale = {
    "Distance": sigma,
    "Energy": epsilon,
    "Mass": m,
    "Time": t,
    "Temperature": T,
    "Pressure": p,
}


def Plot_Inst_Data(
    path_to_data, observable=None, label="Data", cutoff=20000, skip=1, axis=None
):

    if axis is None:
        axis = matplotlyb.pyplot.gca()

    data = numpy.loadtxt(path_to_data)

    if observable is None:
        t = numpy.arange(len(data))
        seaborn.lineplot(x=t[:cutoff:skip], y=data[:cutoff:skip], label=label, ax=axis)
        axis.grid(True)
    else:
        scale = lj_scale[observable]
        timescale = lj_scale["Time"]
        t = numpy.arange(len(data)) * timescale
        seaborn.lineplot(x=t[:cutoff:skip], y=scale * data[:cutoff:skip], label=label, ax=axis)
        axis.grid(True)
