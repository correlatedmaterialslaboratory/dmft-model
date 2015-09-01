from __future__ import division, print_function, unicode_literals


import numpy as np


def density(w_n, E, Sigma, Sigma0):
    """
    Performs the matsubara summation:
         .--           1
       T  >  ---------------------- = N
         '-- iw_n - E - Sigma(iw_n)
    where Sigma(iw_n) asymptotes to Sigma0 at infinity.
    """
    def fermi(x): return 1/(np.exp(x) + 1)

    T = w_n[0]/np.pi            # w_0 = pi * T
    return 2*T * np.sum(1/(1j*w_n - E - Sigma) - 1/(1j*w_n - E - Sigma0)).real + fermi((E + Sigma0)/T)

if __name__ == '__main__':
    pass
