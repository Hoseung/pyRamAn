# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:37:02 2015

@author: hoseung
"""


class Nml():
    def __init__(self, fname=None):
        self.set_fname(fname)
        self.aout = []
        self.x_refine = []
        self.y_refine = []
        self.z_refine = []
        self.r_refine = []

    def loadNml(self, fname=None):
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

