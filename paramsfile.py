from __future__ import division, print_function, unicode_literals


class ParamsFile(object):
    def __init__(self, filename = 'PARAMS'):
        self.filename = filename
        self._fields = [
            ('Sig',            'Sig.out',   str,   "File to output self-energy."),
            ('Gf',             'Gf.out',    str,   "File to output Green's function."),
            ('Delta',          'Delta.inp', str,   "File of input hybridization function."),
            ('cix',            'ctqmc.cix', str,   "Cix file describing atomic states."),
            ('U',              0.0,         float, "Coulomb repulsion (F0)"),
            ('mu',             0.0,         float, "QMC chemical potential"),
            ('beta',           1.0,         float, "Inverse temperature"),
            ('M',              1000000,     int,   "Total number of Monte Carlo steps"),
            ('warmup',         0,           int,   "Warmup number of QMC steps"),
            ('nom',            200,         int,   "Number of Matsubara frequency points sampled"),
            ('Nmax',           2000,        int,   "Maximum perturbation order allowed"),
            ('Segment',        2,           int,   "Whether to use the segment representation of impurity trace"),
            ('tsample',        1000,        int,   "How often to record measurements"),
            ('SampleGtau',     1000,        int,   "How often to update G(tau)"),
            ('GlobalFlip',     1000000000,  int,   "How often to try a global flip"),
            ('CleanUpdate',    1000000,     int,   "How often to make clean update"),
            ('Ncout',          1000000,     int,   "How often to print out info"),
            ('Naver',          1000000000,  int,   "How often to print out debug info"),
            ('Ntau',           1000,        int,   "Number of imaginary time points (only for debugging)"),
            ('aom',            30,          int,   "Number of frequency points used to determine the value of sigma at nom"),
            ('sderiv',         0.02,        float, "Maximum derivative mismatch accepted for tail concatenation"),
            ('FastTrace',      0,           int,   "Speeds up code for full-U calculation"),
            ('LazyTrace',      0,           int,   "Speeds up code for full-U calculation"),
            ('PChangeOrder',   0.9,         float, "Ratio between trial steps: add-remove-a-kink / move-a-kink"),
            ('minM',           1e-10,       float, "The smallest allowed value for the atomic trace"),
            ('minD',           1e-10,       float, "The smallest allowed value for the determinant"),
            ('Ncorrect',       -1,          int,   "Which baths should not be corrected"),
            ('fastFilesystem', 0,           int,   "Whether to write out Gaver and Delta.tau"),
            ]

        #create default attributes
        for (name, default, mytype, doc) in self._fields:
            setattr(self, name, default)

    def readstream(self, stream):
        #remove all white space from a line
        lines = [line.strip() for line in stream.readlines()]
        #ignore lines that are commented out and empty lines
        lines = [line for line in lines if not line.startswith('#') and line != '']

        #output the lines from stream via key and value
        for line in lines:
            print (len(lines))
            record = line.split()
            name = record[0]
            value = record[1]
            mytype = type(getattr(self, name))
            setattr(self, name, mytype(value))
            print (record[0],record[1])

    def read(self, filename = None):
        if not filename:
            filename = self.filename
        with open(filename, 'r') as f:
            self.readstream(f)

    def writestream(self, stream):
        for (name, default, mytype, doc) in self._fields:
            stream.write("{:15s} {:<15} # {}\n".format(name, getattr(self,name), doc))

    def write(self, filename = None):
        if not filename:
            filename = self.filename
        with open(filename, 'w') as f:
            self.writestream(f)
