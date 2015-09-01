from __future__ import division, print_function, unicode_literals


class ParamsFile(object):
    def __init__(self, filename = 'PARAMS'):
        self.filename = filename
        self._fields = [
            ('Sig',          'Sig.out',   str,   "File to output self-energy."),
            ('Gf',           'Gf.out',    str,   "File to output Green's function."),
            ('Delta',        'Delta.inp', str,   "File of input hybridization function."),
            ('cix',          'ctqmc.cix', str,   "Cix file describing atomic states."),
            ('U',            0.0,         float, "Coulomb repulsion (F0)"),
            ('mu',           5.5,         float, "QMC chemical potential"),
            ('beta',         100.0,       float, "Inverse temperature"),
            ('M',            10000000000, int,   "Total number of Monte Carlo steps"),
            ('warmup',       20000000,    int,   "Warmup number of QMC steps"),
            ('tsample',      1000,        int,   "How often to record measurements"),
            ('SampleGtau',   1000,        int,   "How often to update G(tau)"),
            ('GlobalFlip',   1000000000,  int,   "How often to try a global flip"),
            ('CleanUpdate',  1000000,     int,   "How often to make clean update"),
            ('Ncout',        1000000,     int,   "How often to print out info"),
            ('Naver',        1000000000,  int,   "How often to print out debug info"),
            ('Nmax',         2000,        int,   "Maximum perturbation order allowed"),
            ('Ntau',         1000,        int,   "Number of imaginary time points (only for debugging)"),
            ('nom',          200,         int,   "Number of Matsubara frequency points sampled"),
            ('aom',          30,          int,   "Number of frequency points used to determine the value of sigma at nom"),
            ('sderiv',       0.02,        float, "Maximum derivative mismatch accepted for tail concatenation"),
            ('FastTrace',    0,           int,   "Speeds up code for full-U calculation"),
            ('LazyTrace',    0,           int,   "Speeds up code for full-U calculation"),
            ('PChangeOrder', 0.9,         float, "Ratio between trial steps: add-remove-a-kink / move-a-kink"),
            ('minM',         1e-10,       float, "The smallest allowed value for the atomic trace"),
            ('minD',         1e-10,       float, "The smallest allowed value for the determinant"),
            ('Ncorrect',     -1,          int,   "Which baths should not be corrected"),
            ]

        for (name, default, mytype, doc) in self._fields:
            fget = lambda self: getattr(self, '_' + name)
            def fset(self, value):
                assert(type(value) == mytype)
                setattr(self, '_' + name, value)
            setattr(self.__class__, name, property(fget, fset)) # create the property
            setattr(self, '_' + name, default) # create local var to hold value

    def read(self, filename = None):
        if not filename:
            filename = self.filename
        with open(filename, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            lines = [line for line in lines if not line.startswith('#') and line != '']
            print lines
            for line in lines:
                record = line.split()
                print record
                name = record[0]
                value = record[1]
                mytype = type(getattr(self, '_' + name))
                setattr(self, '_' + name, mytype(value))

    def write(self, filename = None):
        if not filename:
            filename = self.filename
        with open(filename, 'w') as f:
            for (name, default, mytype, doc) in self._fields:
                f.write("{:15s} {:<15} # {}\n".format(name, getattr(self, '_' + name), doc))
