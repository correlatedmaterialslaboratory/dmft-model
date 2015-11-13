#!/usr/bin/env python


from __future__ import division, print_function


import pythtb
import pytriqs.gf.local as gf
import ConfigParser


def parse_configs(filename):
    config = ConfigParser.ConfigParser()
    config.read(filename)
    return config
    # How do I want to access parameters? Just leave them inside the config parser?

def read_self_energy(Sinfty):
    indices = ['1', '2', '3', '4']
    Sigma = gf.GfImFreq(indices = indices, beta = 100, name = '$\Sigma_{ij}^\mathrm{imp}$')
    for i in indices:
        Sigma[i,i] = Sinfty
    return Sigma

def main():
    params = parse_configs('dmft.cfg')
    U = params.getfloat('Physical Parameters', 'U')
    nd0 = params.getfloat('Physical Parameters', 'nd0')
    Edc = U*(nd0 - 0.5)
    # H = construct_model()
    Sigma = read_self_energy(Edc)
    print(Sigma.data[0])

if __name__ == '__main__':
    main()
