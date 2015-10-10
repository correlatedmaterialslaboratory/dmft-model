from __future__ import division, print_function, unicode_literals


import unittest
import io
from os import path
from paramsfile import ParamsFile


class ParamsFileTest(unittest.TestCase):
    def setUp(self):
        CWD = path.dirname(__file__)
        self.fparams = path.join(CWD, 'PARAMS')
        self.fparams_blank = path.join(CWD, 'PARAMS.blank')

    def assertValues(self, pf):
        expected = [
            ('Sig',             'Sig.out'),
            ('Gf',              'Gf.out'),
            ('Delta',           'Delta.inp'),
            ('cix',             '12actqmc.cix'),
            ('U',               0.0),
            ('mu',              7.8),
            ('beta',            100.0),
            ('M',               200000000),
            ('warmup',          30000000),
            ('tsample',         1000),
            ('SampleGtau',      1000),
            ('GlobalFlip',      1000000000),
            ('CleanUpdate',     1000000),
            ('Ncout',           1000000),
            ('Naver',           1000000000),
            ('Nmax',            2000),
            ('Ntau',            1000),
            ('nom',             100),
            ('aom',             7),
            ('sderiv',          0.02),
            ('FastTrace',       0),
            ('LazyTrace',       0),
            ('PChangeOrder',    0.9),
            ('minM',            1e-10),
            ('minD',            1e-10),
            ('Ncorrect',        -1),
            ('fastFilesystem',  1),
            ('Segment',         2),
            ]

        for attr,val in expected:
            self.assertEqual(getattr(pf, attr), val)

    def assertValuesDefaults(self, pf):
        defaults = [
            ('Sig',            'Sig.out'),
            ('Gf',             'Gf.out'),
            ('Delta',          'Delta.inp'),
            ('cix',            'ctqmc.cix'),
            ('U',              0.0),
            ('mu',             0.0),
            ('beta',           1.0),
            ('M',              1000000),
            ('warmup',         0),
            ('nom',            200),
            ('Nmax',           2000),
            ('Segment',        2),
            ('tsample',        1000),
            ('SampleGtau',     1000),
            ('GlobalFlip',     1000000000),
            ('CleanUpdate',    1000000),
            ('Ncout',          1000000),
            ('Naver',          1000000000),
            ('Ntau',           1000),
            ('aom',            30),
            ('sderiv',         0.02),
            ('FastTrace',      0),
            ('LazyTrace',      0),
            ('PChangeOrder',   0.9),
            ('minM',           1e-10),
            ('minD',           1e-10),
            ('Ncorrect',       -1),
            ('fastFilesystem', 0),
            ]

        for attr,val in defaults:
            self.assertEqual(getattr(pf, attr), val)

    def test_read(self):
        pf = ParamsFile()
        pf.read(self.fparams)
        self.assertValues(pf)

    def test_readblank(self):
        pf = ParamsFile(self.fparams_blank)
        pf.read()
        self.assertValuesDefaults(pf)

    def test_write(self):
        pf0 = ParamsFile(self.fparams)
        pf0.read()
        out = io.StringIO()
        pf0.writestream(out)
        out.seek(0)
        pf1 = ParamsFile()
        pf1.readstream(out)
        self.assertValues(pf1)

    def test_writeblank(self):
        pf0 = ParamsFile(self.fparams_blank)
        pf0.read()
        out = io.StringIO()
        pf0.writestream(out)
        out.seek(0)
        pf1 = ParamsFile()
        pf1.readstream(out)
        self.assertValuesDefaults(pf1)
