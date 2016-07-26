# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:10:48 2015

@author: hoseung
"""
import numpy as np

class Info:
    def __init__(self, nout = None, base = './',
                 fn = None, load=True, data_dir=None,
                 cosmo=True):
        if data_dir is None:
            from general import defaults
            df = defaults.Default()
            data_dir = df.dir_snapshot
        self.data_dir = data_dir
        self.fn = fn
        self.cosmo=cosmo
        
        if nout is not None or base is not None or fn is not None:
            self.setup(nout=nout, base=base, fn=fn)
        if load:
            self.read_info()
        
    def setup(self, nout = None, base = './', fn = None):
        try:  # set nout
            self._set_nout(nout)
        except:
            print("info: NOUT is not given")

        try: # set base directory
            self._set_base(base)
        except:
            print("info: BASE is not given")

        self.update_fn(fn)
 
        if self.all_set():
            self.read_info()  # sets ncpu_tot,


    def __call__(self, *args):
        # Function emulation      
        return self.__init__(*args)

    def all_set(self):
        return (self.nout is not None) & (self.base is not None)

    def update_fn(self, fn=None):
        import os
        try:
            snout = str(self.nout).zfill(5)
            self.fn = os.path.join(self.base, self.data_dir) + 'output_' + snout + '/info_' + snout + '.txt'
        except:
            if self.nout is None and self.base is None:
                try:
                    self.set_fn(fn)
                except:
                    print("Info file name is not given. ")
                    print("Current value is :", self.fn)

    def _set_base(self, base):
        self.base = base
        
    def _set_nout(self, nout):
        self.nout = nout
        self.snout = str(self.nout).zfill(5)

    def _set_ranges(self, ranges):
        self.ranges = ranges

    def _set_cpus(self, cpus):
        self.cpus = np.array(cpus)

    def set_fn(self, fn, verbose=False):
        self.fn = fn
        if verbose:
            print("Info file name: %s" % self.fn)

# Weak "internal use" indicator.
# from M import * does not import _blah_blah.

    def _cal_units(self, arr, rarr):
        import utils.cosmology
        # in cgs unit
        kpc = 3.08e21
        twopi = 6.2831853e0
        hplanck = 6.6262000e-27
        eV = 1.6022000e-12
        kB = 1.38e-16
        clight = 2.9979250e+10
        Gyr = 3.1536000e+16
        X = 0.76
        Y = 0.24
        rhoc = 1.8800000e-29
        mH = 1.6600000e-24
        mu_mol = 1.2195e0
        G = 6.67e-8
        m_sun = 1.98892e33
        scale_l = rarr[8]
        scale_d = rarr[9]
        scale_t = rarr[10]
        scale_v = scale_l / scale_t
        scale_T2 = mH/kB * scale_v**2
        scale_nH = X/mH * scale_d
        scale_Z = 1./0.02
        scale_flux = scale_v * scale_d * kpc * kpc * Gyr/m_sun

        self.ncpu_tot = arr[0]
        self.boxtokpc = rarr[0] * scale_l/kpc
        self.boxlen = rarr[0]
        self.lmin = arr[2]
        self.lmax = arr[3]
        self.unit_l = scale_l
        self.unit_d = scale_d
        self.unit_t = scale_t
        self.unit_v = scale_v
        self.unit_nH = scale_nH
        self.unit_T2 = scale_T2
        self.unit_Z = scale_Z
        self.kms = scale_v/1e5
        self.unit_flux = scale_d * scale_v * (1e-9*Gyr)*kpc/m_sun
        self.punit_m = scale_d * scale_l**3/m_sun
        self.pboxsize = rarr[0] * scale_l/(kpc*1000.)
        self.time = rarr[1]
        self.aexp = rarr[2]
        self.zred = 1/rarr[2]-1
        self.H0 = rarr[3]
        self.h = self.H0 * self.aexp
        self.om = rarr[4]
        self.ol = rarr[5]
        self.ok = rarr[6]
        self.ob = rarr[7]
        self.msun = scale_d*scale_l**3/m_sun
        self.cboxsize = self.H0 * self.pboxsize / self.aexp * 1e-2
        if self.cosmo:
            self.tGyr = utils.cosmology.time2gyr(rarr[1], z_now = self.zred, info = self)

    def keys(self):
        from pprint import pprint
        pprint(vars(self))

    def read_info(self, base=None, nout=None, verbose=False):
        """ backward compatibility. but use .load instead"""
        self.load(base=base, nout=nout, verbose=verbose)
        
    def load(self, base=None, nout=None, verbose=False):
        """
            parameters: nout, base
        """
#        print("Loading INFO")
        if self.base is None:
            if base is None:
                raise ValueError("Working directory is not determined")
            else:
                self._set_base(base)

        if self.nout is None:
            if nout is None:
                raise ValueError("nout is not determined")
            else:
                self._set_nout(nout)

        if self.fn is None:
            self.update_fn()
        if verbose:
            print(self.fn)

        with open(self.fn) as f:
            arr1 = []
            arr2 = []

            # Parse info file. 
            for i in range(5):
                arr1.append(int(str.split(f.readline(), '=')[1]))
            self.ndim = arr1[1]
            self.nstep_coarse = int(str.split(f.readline(), '=')[1])
            f.readline()  # empty line
            for i in range(11):
                arr2.append(float(str.split(f.readline(), '=')[1]))
            self._cal_units(arr1, arr2)
            self.hilbertkey = np.zeros((2, self.ncpu_tot + 1),
                                       dtype=np.float64)
            f.readline()
            f.readline()
            f.readline()
            for i in range(self.ncpu_tot):
                keystr = (str.split(f.readline()))[1:]
                self.hilbertkey[:, i] = keystr[0]
                self.hilbertkey[:, i + 1] = keystr[1]


class Nml():
    def __init__(self, fname=None, setup=False):
        if fname is not None:
            # If filename is given, laod it by default
            setup = True
            self.set_fname(fname)
        self.aout = []
        self.x_refine = []
        self.y_refine = []
        self.z_refine = []
        self.r_refine = []
        if setup:
            self.load(fname=self.fname)

    def load(self, fname=None):
        import numpy as np
        import re
        if fname is None:
            fname = self.fname
        with open(fname, mode='r') as f:
            for line in f.readlines():
                if 'aout' in line:
                    self.aout = np.asarray([float(i) for i in re.split("=|,",line)[1:]])
                if 'x_refine' in line:
                    tmp = np.asarray([float(i) for i in re.split("=|,|\*",line)[1:]])
                    self.x_refine = tmp[tmp < 1]
                if 'y_refine' in line:
                    tmp = np.asarray([float(i) for i in re.split("=|,|\*",line)[1:]])
                    self.y_refine = tmp[tmp < 1]
                if 'z_refine' in line:
                    tmp = np.asarray([float(i) for i in re.split("=|,|\*",line)[1:]])
                    self.z_refine = tmp[tmp < 1]
                if 'r_refine' in line:
                    tmp = np.asarray([float(i) for i in re.split("=|,|\*",line)[1:]])
                    self.r_refine = tmp[tmp < 1]

    def set_fname(self, fname):
        self.fname = fname


class RefineParam():
    def __init__(self, name=None):
        self.name = name
        self.nnout = 0
        self.nzoomlevel = 0
        self.aexp = []
        self.x_refine = []
        self.y_refine = []
        self.z_refine = []
        self.r_refine = []

    def loadRegion(self, fname, aexp=None):
        """
        """
        with open(fname, mode='r') as f:
            f.readline()  # nnout
            self.nnout = int(f.readline())
            f.readline()  # nzoomlevel
            self.nzoomlevel = f.readline()
            f.readline()  # AEXP
            self.aexp = [float(i) for i in f.readline().split(",")]
            f.readline()  # X_REFINE
            for i in range(self.nnout):
                self.x_refine.append(float(f.readline()))
            f.readline()  # Y_REFINE
            for i in range(self.nnout):
                self.y_refine.append(float(f.readline()))
            f.readline()  # Z_REFINE
            for i in range(self.nnout):
                self.z_refine.append(float(f.readline()))
            f.readline()  # R_REFINE
            for i in range(self.nnout):
                self.r_refine.append(float(f.readline()))

