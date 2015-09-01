#!/usr/bin/env python

import optparse
import scipy as sp
from scipy import interpolate

def main():
    usage = """usage: %prog [options]

    Interpolates existing self-energy on Matsubara mesh for a new temperature.
    """

    parser = optparse.OptionParser(usage)

    parser.add_option("-i", "--insig", default = 'sig.inp', help = "File containing original self-energy.")
    parser.add_option("-o", "--outsig", default = 'sig.out', help = "File to write interpolated self-energy.")
    parser.add_option("-b", "--beta", type = "float", default = 100, help = "Inverse temperature (beta).")

    (opts, args) = parser.parse_args()

    print 'insig = %s, outsig = %s, beta = %s' %  (opts.insig, opts.outsig, opts.beta)

    # Read the rest of the input file
    sig_old = sp.loadtxt(opts.insig).T
    w_old = sig_old[0]
    w_max = w_old[-1]
    print 'omega(old) = ', w_old

    # Create new Matsubara mesh
    nmax = (opts.beta * w_max / sp.pi - 1) / 2.
    w_new = sp.pi * (2*sp.arange(nmax) + 1) / opts.beta
    print 'omega(new) = ', w_new

    # Interpolate each component of old self-energy onto new mesh
    sig_new=[]
    sig_new.append(w_new)
    for i in range(1,len(sig_old)):        
        tck = interpolate.splrep(w_old, sig_old[i], k=3, s=0)
        newy = interpolate.splev(w_new, tck)
        sig_new.append(newy)
    sig_new = sp.array(sig_new)

    print 'Output written to: ', opts.outsig
    sp.savetxt(opts.outsig, sig_new.T)

if __name__=='__main__':
    main()
