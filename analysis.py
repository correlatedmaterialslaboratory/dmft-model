from __future__ import division, print_function, unicode_literals


import numpy as np


def _fit_three_points(w_n, ImSigma):
    x = w_n[:3]
    y = ImSigma[:3]
    z = np.polyfit(x, y, 2)
    p = np.poly1d(z)
    return p

def extract_z(w_n, ImSigma):
    p = _fit_three_points(w_n, ImSigma)
    pderiv = np.polyder(p)
    z = 1/(1-pderiv(0))
    return z

def extract_scattering(w_n, ImSigma):
    '''
    Scattering rate is -Im Sigma(w = 0)
    '''
    p = _fit_three_points(w_n, ImSigma)
    return -p(0)
