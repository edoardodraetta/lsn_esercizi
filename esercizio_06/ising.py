import numpy as np


def energy(J, Ns, start, stop, N_Pts):
    ''' 1D ising spin model energy '''
    T = np.linspace(start, stop, N_Pts)
    th = np.tanh(J / T)
    thN = th ** Ns
    ch = 1 / th
    e = -J * (th + ch * thN) / (1 + thN)

    return e


def heat(J, Ns, start, stop, N_Pts):
    T = np.linspace(start, stop, N_Pts)
    beta = 1 / T
    th = np.tanh(J / T)
    thN = th ** Ns
    ch = 1 / th
    heat = ((beta * J) ** 2) * (
        ((1 + thN + (Ns - 1) * (th ** 2) + (Ns - 1) *
          (ch ** 2) * thN) / (1 + thN)
         ) - Ns * ((th + ch * thN) / (1 + thN)) ** 2)

    return heat


def susceptibility(J, Ns, start, stop, N_Pts):
    T = np.linspace(start, stop, N_Pts)
    beta = 1 / T
    th = np.tanh(J / T)
    thN = th ** Ns
    X = beta * np.exp(2 * beta * J) * (1 - thN) / (1 + thN)

    return X


def magnetization(h, J, Ns, start, stop, N_Pts):
    T = np.linspace(start, stop, N_Pts)
    h = 0.02  # external field
    b = 1 / T
    l1 = np.exp(b*J)*np.cosh(b*h)+np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J))
    l2 = np.exp(b*J)*np.cosh(b*h)-np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J))
    Z = l1**Ns + l2**Ns
    M = (np.exp(b*J)*np.sinh(b*h)*((l1**(Ns-1))*(1+np.exp(b*J)*np.cosh(b*h)/np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J)))
            + (l2**(Ns-1))*(1-np.exp(b*J)*np.cosh(b*h)/np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J)))))/(Z)

    return M
