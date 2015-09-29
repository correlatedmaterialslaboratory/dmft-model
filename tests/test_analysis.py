from __future__ import division, print_function, unicode_literals


import unittest

from math import pi
import numpy as np
from analysis import extract_z, extract_scattering


class SigmaTest(unittest.TestCase):
    def test_quad(self):
        T = 0.01
        w_n = pi * (2*np.arange(10)+1) * T
        ImSigma = w_n * (w_n - 10)
        z = extract_z(w_n, ImSigma)
        self.assertAlmostEqual(z, 1/11)
        gamma = extract_scattering(w_n, ImSigma)
        self.assertAlmostEqual(gamma, 0)

    def test_linear(self):
        T = 0.02
        w_n = pi * (2*np.arange(10)+1) * T
        ImSigma = -0.5 - w_n
        z = extract_z(w_n, ImSigma)
        self.assertAlmostEqual(z, 1/2)
        gamma = extract_scattering(w_n, ImSigma)
        self.assertAlmostEqual(gamma, 0.5)
