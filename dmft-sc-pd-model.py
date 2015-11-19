#!/usr/bin/env python


import pythtb as ptb            # Sinisa's and David's Python tight-binding module
import numpy as np
import os, sys, shutil, subprocess, glob
import argparse
import ConfigParser

from numpy import pi

class pd_model(object):
    def __init__(self, Ed, Ep, tpd, tpp, tppp):
        self.Ed = Ed
        self.Ep = Ep
        self.tpd = tpd
        self.tpp = tpp
        self.tppp = tppp

        dim_k = 2
        dim_r = 2
        lat = [(1,0), (0,1)]
        orb = [(0,0), (.5,0), (0,.5)]

        self.icorr_orbs = [0]   # list of indices of correlated orbitals

        # composition: has-a relationship with tb_model
        self.model = ptb.tb_model(dim_k, dim_r, lat, orb)
        self.model.set_onsite([Ed,Ep,Ep])
        
        # (amplitude, i, j, [latt. vector to cell containing j])
        hoppings = [
            ( tpd, 0, 1, (0,0)),     # p-d hoppings
            (-tpd, 0, 2, (0,0)),
            (-tpd, 0, 1, (-1,0)),
            ( tpd, 0, 2, (0,-1)),
            (-tpp, 1, 2, (0,0)),     # diagonal p-p hoppings
            (-tpp, 1, 2, (1,-1)),
            ( tpp, 1, 2, (1,0)),
            ( tpp, 1, 2, (0,-1)),
            ( tppp,1, 1, (1,0)),     # tpp'
            ( tppp,2, 2, (0,1)),
            ]

        for amp,i,j,lat in hoppings:
            self.model.set_hop(amp, i, j, lat)

    @property
    def Norbs(self):
        return self.model.get_num_orbitals()

    @property
    def Ncorr(self):
        return len(self.icorr_orbs)

    @property
    def corr_orbs(self):
        # compute indices of correlated orbitals
        orbs = self.model.get_orbitals() # list of orbital coordinates
        return [orbs[i] for i in self.icorr_orbs]

    def make_supercell(self, *args, **kwargs):
        sc_model = self.model.make_supercell(*args, **kwargs) # 2x2 supercell
        # create new pd_model with the sc_model inside
        new = pd_model(self.Ed, self.Ep, self.tpd, self.tpp, self.tppp)
        new.model = sc_model
        [sc_red_lat] = args      # matrix M[i,j]; i indexes the new vectors, j indexes the coordinates
        M = np.array(sc_red_lat).T
        reduced_coords = [[x % 1.0 for  x in np.dot(M, orb)] for orb in sc_model.get_orbitals()]
        d_orbital = (0.0,0.0)
        new.icorr_orbs = [i for (i,coord) in enumerate(reduced_coords) if tuple(coord) == d_orbital]
        return new

    def dos(self, Nkpts = 100, rand = False):
        # k-points for first-quadrant only
        if rand:
            kpts = np.random.rand(Nkpts**2, 2) * 0.5 # sample BZ randomly
        else:
            kpts = [((i+.5)/(2*Nkpts), (j+.5)/(2*Nkpts)) for i in range(Nkpts) for j in range(Nkpts)]
        bands = self.model.solve_all(kpts)

        # histogram points intelligently
        # |   |   |   | ... |   |
        #   ^                 ^
        #   |                 |
        #  Emin              Emax
        Emin = np.min(bands)
        Emax = np.max(bands)
        Nbins = int((Nkpts**2)/33.)  # 33 points is a statistically "large" sample
        dE = (Emax - Emin)/(Nbins-1)
        # need one extra bin because dos[right_bin] will otherwise fail for Emax
        dos = np.zeros(Nbins+1)
    
        for E in bands.flat:
            ibin = (E-Emin)/dE
            left_bin = int(ibin)
            right_bin = left_bin + 1
            right_weight = ibin % 1
            left_weight = 1-right_weight
            dos[left_bin] += left_weight
            dos[right_bin] += right_weight

        Emesh = np.linspace(Emin, Emax, Nbins)
    
        return Emesh, dos[:-1]/(Nkpts**2)

    def display(self):
        print "Number of orbitals =", self.Norbs
        print "Number of correlated orbitals =", self.Ncorr
        print "Coordinates of correlated orbitals =", self.corr_orbs
        print "Indices of correlated orbitals =", self.icorr_orbs
        print self.model.display()


def mat2str(m, decimals = 5):
    format_one = '{:' + str(decimals + 3) + '.' + str(decimals) + 'f}'
    format = format_one + ' ' + format_one
    return "\n".join(["  ".join([format.format(e.real, e.imag) for e in row]) for row in m])


def read_spectral(filename, matrix_form):
    data = np.loadtxt(filename).T
    omega_n = data[0]
    components = data[1::2] + 1j * data[2::2]
    dim_i, dim_j = np.shape(matrix_form)
    ret = np.zeros((len(omega_n), dim_i, dim_j), dtype = complex)
    print ret.shape
    for i in range(dim_i):
        for j in range(dim_j):
            index = matrix_form[i,j]
            if index > 0:
                ret[:,i,j] = components[index-1]
    return omega_n, ret


def write_spectral(filename, omega_n, data):
    data_flat = np.empty((len(omega_n), 6)) # 3 unique components * 2 for real/imag
    data_flat[:,0::2] = data[:,0,[0,1,3]].real
    data_flat[:,1::2] = data[:,0,[0,1,3]].imag
    np.savetxt(filename, np.hstack((omega_n[:,np.newaxis], data_flat)))


def read_self_energy(filename, matrix_form, Nw, T, Sinfty):
    if os.path.isfile(filename):
        omega_n, Sigma = read_spectral(filename, matrix_form)
    else:
        dim_i, dim_j = matrix_form
        Sigma = np.zeros((Nw, dim_i, dim_j), dtype = complex)
        assert(dim_i == dim_j)
        for i in range(dim_i):
            Sigma[:,i,i] = Sinfty
        omega_n = pi * (2*np.arange(Nw)+1) * T # w_n = pi (2n+1) T
        write_spectral(filename, omega_n, Sigma)

    print 'omega_n =', omega_n[:10], "..."
    print "shape(Sigma) =", Sigma.shape
    print "Sigma(pi*T) ="
    print mat2str(Sigma[0])

    return omega_n, Sigma


def embed_self_energy(pd, Sigma):
    # embed self-energy into "lattice"
    Nw, _ ,_ = Sigma.shape
    Norbs = len(pd.model.get_orbitals())
    Sigma_c = np.zeros((Nw, Norbs, Norbs), dtype = complex)
    for i,ic in enumerate(pd.icorr_orbs):
        for j,jc in enumerate(pd.icorr_orbs):
            Sigma_c[:,ic,jc] = Sigma[:,i,j]

    print "shape(Sigma_c) =", Sigma_c.shape
    print "Sigma_c(pi*T) ="
    print mat2str(Sigma_c[0])

    return Sigma_c


def construct_klist(Nk):
    kmesh = np.linspace(-.5+.5/Nk, .5-.5/Nk, Nk)
    return [(kx,ky) for kx in kmesh for ky in kmesh]


def compute_phases(pd):
    # We choose the tight-binding convention where orbitals have fractional
    # coordinates within the unit cell (as opposed to omitting the internal
    # coordinates). This implies phase factors when perform sums over the sites
    # within the unit cell.
    #
    # tdiff[i,j] = i * k . (\tau_i - \tau_j))
    #   \tau_i is coordinate of ith orbital, and k in [-pi,pi] is BZ momentum
    tt = np.array(pd.corr_orbs).T
    ti = tt[:,:,np.newaxis]
    tj = tt[:,np.newaxis,:]
    tdiff = 1j * 2*pi * (ti - tj) # need explicit 2pi factor since k-points are in [-.5,.5]
    print "orbital locations ="
    print tt
    print "shape(tdiff) =", tdiff.shape
    print "i(t_i-t_j) ="
    print mat2str(tdiff[0])
    print
    print mat2str(tdiff[1])
    return tdiff


def compute_local_greens_function(pd, omega_n, Sigma_c, Nk):
    # Compute cluster Green's function from "lattice"
    #   explicitly loop over k-points
    #   {i,j} and omega_n handled by broadcasting
    Nw, Norbs, Norbs = Sigma_c.shape
    Ncorr = pd.Ncorr
    kpts = construct_klist(Nk)

    tdiff = compute_phases(pd)
    Gc = np.zeros((Nw, Ncorr, Ncorr), dtype = complex)
    print "shape(Gc) =", Gc.shape
    I = np.identity(Norbs)[np.newaxis,:,:]
    iw = 1j * omega_n[:,np.newaxis,np.newaxis] * I # broadcasted Matsubara frequencies
    for ik,k in enumerate(kpts):
        print ik, k
        Hk = pd.model.ham(k)[np.newaxis,:,:]
        G_bare = np.array([np.linalg.inv(Ginv)[pd.icorr_orbs,:][:,pd.icorr_orbs] for Ginv in iw - Hk - Sigma_c])
        phases = np.exp(np.tensordot(k, tdiff, 1))[np.newaxis,:,:]
        Gc += G_bare * phases   # element-wise multiplication
    Gc /= len(kpts)
    print "Gc[0] ="
    print mat2str(Gc[0])
    write_spectral("Gc.inp", omega_n, Gc)

    return Gc


def compute_greens_function(pd, omega_n, Sigma_c, Nk):
    # Compute full cluster Green's function from "lattice"
    #   explicitly loop over k-points
    #   {i,j} and omega_n handled by broadcasting
    Nw, Norbs, Norbs = Sigma_c.shape
    Ncorr = pd.Ncorr
    kpts = construct_klist(Nk)

    tt = np.array(pd.model.get_orbitals()).T
    ti = tt[:,:,np.newaxis]
    tj = tt[:,np.newaxis,:]
    tdiff = 1j * 2*pi * (ti - tj) # need explicit 2pi factor since k-points are in [-.5,.5]

    Gc = np.zeros((Nw, Norbs, Norbs), dtype = complex)
    print "shape(Gc) =", Gc.shape
    I = np.identity(Norbs)[np.newaxis,:,:]
    iw = 1j * omega_n[:,np.newaxis,np.newaxis] * I # broadcasted Matsubara frequencies
    for ik,k in enumerate(kpts):
        print ik, k
        Hk = pd.model.ham(k)[np.newaxis,:,:]
        G_bare = np.array([np.linalg.inv(Ginv) for Ginv in iw - Hk - Sigma_c])
        phases = np.exp(np.tensordot(k, tdiff, 1))[np.newaxis,:,:]
        Gc += G_bare * phases   # element-wise multiplication
    Gc /= len(kpts)
    print "Gc[0] ="
    print mat2str(Gc[0])
    # write_spectral("Gc.inp", omega_n, Gc)

    for n in range(Nobs):
        matsubara.density(Gc[:,n,n])

    return Gc


def compute_impurity_levels(pd, Nk):
    # compute impurity levels
    Eimp = np.zeros((pd.Ncorr, pd.Ncorr), dtype = complex)
    kpts = construct_klist(Nk)
    tdiff = compute_phases(pd)
    for k in kpts:
        Hk = pd.model.ham(k)[pd.icorr_orbs,:][:,pd.icorr_orbs]
        phases = np.exp(np.tensordot(k, tdiff, 1))
        Eimp += Hk * phases     # element-wise multiplication
    Eimp /= len(kpts)
    print "Eimp ="
    print mat2str(Eimp)
    return Eimp


def compute_hybridization(omega_n, Gc, Eimp, Sigma):
    # compute hybridization function
    Ncorr = len(Eimp)
    I = np.identity(Ncorr)[np.newaxis,:,:]
    iw = 1j * omega_n[:,np.newaxis,np.newaxis] * I # broadcasted Matsubara frequencies
    Ginv = np.array([np.linalg.inv(G) for G in Gc])
    Delta = iw - Eimp[np.newaxis,:,:] - Ginv - Sigma
    print "shape(Delta) =", Delta.shape
    print "Delta(pi*T) ="
    print mat2str(Delta[0])
    return Delta
    

def mix_delta(iter, Delta, delta_mixing = 0.0):
    cur_Delta = Delta
    if iter > 0:
        _, prev_Delta = read_spectral("Delta."+str(iter-1))
    else:
        prev_Delta = Delta
    return cur_Delta * (1 - delta_mixing) + prev_Delta * delta_mixing

def save_spectral_quantities(iter, omega_n, Gc, Delta, Eimp):
    write_spectral("Gc."+str(iter), omega_n, Gc)
    write_spectral("Delta."+str(iter), omega_n, Delta)
    f = open("Eimp."+str(iter), 'w')
    f.write(mat2str(Eimp))
    f.close()


def solve_impurity(omega_n, Eimp, Delta, remove_status_files = False):
    # user must have setup PARAMS, cix file and mpi_prefix
    write_spectral("Delta.inp", omega_n, Delta)

    # delete old status files
    if remove_status_files:
        for filename in glob.glob('status.*'):
            os.remove(filename)

    # execute CTQMC
    mpi_prefix = open("mpi_prefix", 'r').read().strip().split()
    cmd = mpi_prefix + "ctqmc".split()
    print ' '.join(cmd)
    with open('log.ctqmc', 'w') as logfile:
        ret = subprocess.call(cmd, stdout = logfile, stderr = logfile)
        if ret != 0:
            print "Error in CTQMC. Check log.ctqmc for error message."

    # gaussian smooth the self-energy to remove noise at intermediate-frequencies
    sigma_filename = 'Sig.out'
    cmd = ['broad', '-w', '0.0', '-k', '0.2', sigma_filename]
    with open(sigma_filename + '.smoothed', 'w') as outfile:
        ret = subprocess.call(cmd, stdout = outfile, stderr = sys.stdout)
        if ret != 0:
            print "Error broadening self-energy."

    omega_n, Sigma = read_spectral(sigma_filename + '.smoothed')
    return Sigma


def save_impurity_output(iter, omega_n, Sigma):
    suffix = '.' + str(iter)
    write_spectral("Sigma" + suffix, omega_n, Sigma)
    to_copy = ['log.ctqmc', 'histogram.dat', 'Probability.dat', 'Gf.out', 'Sig.out']
    for filename in to_copy:
        shutil.copy(filename, filename + suffix)

def parse_commandline():
    parser = argparse.ArgumentParser(description = 'cluster-DMFT for superconducting p-d model')
    parser.add_argument('-c', dest = 'config_file', type=str, default='dmft.cfg', help = 'Configuration file')
    return parser.parse_args()


class ConfigFile:
    def __init__(self):
        self._config = ConfigParser.ConfigParser()
        self._fields = [
            ('T',         1.0,  float, "Physical Parameters",  "Temperature (eV)"),
            ('mu',        0.0,  float, "Physical Parameters",  "Chemical potential (eV)"),
            ('Ed',        0.0,  float, "Physical Parameters",  "Onsite d-orbital energy (eV)"),
            ('Ep',        0.0,  float, "Physical Parameters",  "Onsite p-orbital energy (eV)"),
            ('tpd',       0.0,  float, "Physical Parameters",  "p-d hopping amplitude (eV)"),
            ('tpp',       0.0,  float, "Physical Parameters",  "p-p hopping amplitude (eV)"),
            ('tppp',      0.0,  float, "Physical Parameters",  "p-p hopping amplitude (eV)"),
            ('U',         0.0,  float, "Physical Parameters",  "Onsite Coulomb repulsion (eV)"),
            ('nd0',       0.0,  float, "Physical Parameters",  "Nominal valence for double-counting"),
            ('Nw',        1000, int,   "Numerical Parameters", ""),
            ('Nk',        20,   int,   "Numerical Parameters", ""),
            ('max_iter',  100,  int,   "Numerical Parameters", ""),
            ('delta_mix', 0.75, float, "Numerical Parameters", "How much of previous Delta to mix into new"),
            ('supercell', [[1,0], [0,1]], self._intmatrix, "DMFT", "Basis vectors for supercell"),
            ('test',      [[0]], self._intmatrix, "DMFT", "DOC"),
            ]
        # create defaults
        for (name, default, typefunc, section, doc) in self._fields:
            setattr(self, name, default)

    def _intmatrix(self, string):
        matrix = []
        for row in string.split('\n'):
            if row.strip():     # if row is not empty whitespace
                matrix.append([int(x) for x in row.split()])
        return matrix

    def read(self, filename = 'dmft.cfg'):
        self._config.read(filename)
        for (name, default, typefunc, section, doc) in self._fields:
            if self._config.has_option(section, name):
                setattr(self, name, typefunc(self._config.get(section, name)))

    def write(self):
        pass

def main():
    args = parse_commandline()
    cfg = ConfigFile()
    cfg.read(args.config_file)

    # chemical potential incorporated into onsite energies of model
    # double counting goes into d-level
    Edc = cfg.U*(cfg.nd0 - 0.5)         # fixed double-counting
    pd = pd_model(cfg.Ed-cfg.mu-Edc, cfg.Ep-cfg.mu, cfg.tpd, cfg.tpp, cfg.tppp)
    pd2x2 = pd.make_supercell(cfg.supercell)
    pd2x2.display()

    iter = 0
    omega_n, Sigma = read_self_energy("Sigma.inp", pd2x2, cfg.Nw, cfg.T, Edc)
    while iter < cfg.max_iter and not os.path.isfile('STOP'):
        Sigma_c = embed_self_energy(pd2x2, Sigma)
        Gc = compute_local_greens_function(pd2x2, omega_n, Sigma_c, Nk)
        Eimp = compute_impurity_levels(pd2x2, Nk)
        Delta = compute_hybridization(omega_n, Gc, Eimp, Sigma)
        Delta = mix_delta(iter, Delta, delta_mixing)
        save_spectral_quantities(iter, omega_n, Gc, Delta, Eimp)
        Sigma = solve_impurity(omega_n, Eimp, Delta, remove_status_files = False)
        iter += 1
        save_impurity_output(iter, omega_n, Sigma)

if __name__ == '__main__':
    main()
