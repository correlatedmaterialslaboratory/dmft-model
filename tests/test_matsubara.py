from __future__ import division, print_function, unicode_literals


import unittest

import numpy as np
from matsubara import density


class DensityTest(unittest.TestCase):
    def matsubara(self, T, N = 1000):
        return np.pi * (2*np.arange(N) + 1) * T 

    def fermi(self, x):
        return 1/(np.exp(x) + 1)

    def test_zero(self):
        T = 0.01
        w_n = self.matsubara(T)
        E = 0
        Sigma = w_n * 0
        Sigma0 = 0
        N = density(w_n, E, Sigma, Sigma0)
        self.assertEqual(N, 0.5)

    def test_pole(self):
        T = 0.1
        w_n = self.matsubara(T)
        E = 0.3
        Sigma = w_n * 0
        Sigma0 = 0
        N = density(w_n, E, Sigma, Sigma0)
        self.assertEqual(N, self.fermi(E/T))

    def test_sigma_pole(self):
        T = 0.05
        w_n = self.matsubara(T, N = 10000)
        E = 0.3
        V2 = 0.5
        Sigma0 = 0.314159
        Sigma = V2/(1j*w_n) + Sigma0
        N = density(w_n, E, Sigma, Sigma0)

        a = (E+Sigma0 + np.sqrt((E+Sigma0)**2 + 4*V2))/2
        b = (E+Sigma0 - np.sqrt((E+Sigma0)**2 + 4*V2))/2
        Nref = a/(a-b) * self.fermi(a/T) - b/(a-b) * self.fermi(b/T)

        self.assertAlmostEqual(N, Nref)


if __name__ == "__main__":
    unittest.main()
