# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 20:40:35 2015

@author: hoseung
"""

class Hmgal():
    def __init__(self, fn):
        self.nbodies=0
        self.massp=0.0
        self.aexp=0.0
        self.omegat=0
        self.age=0.0
        self.halnum=0
        self.subnum=0
        self.fn = fn
        self.f = open(fn, "rb")
        self.load_meta()
        self.load_data()
        self.f.close()
        
        def load_meta(self):
            import numpy as np
            from load.utils import read_fortran
            f = self.f
            self.nbodies = read_fortran(f, np.dtype('i4'), 1)
            self.massp = read_fortran(f, np.dtype('f4'), 1)
            self.aexp = read_fortran(f, np.dtype('f4'), 1)
            self.omegat = read_fortran(f, np.dtype('f4'), 1)
            self.age = read_fortran(f, np.dtype('f4'), 1)
            self.ngal, self.nsubgal = read_fortran(f, np.dtype('i4'), 2)

        def load_data(self, tothal=None):
            f = self.f
            import numpy as np
            from load.utils import read_fortran
            dtype_gal = [('np', '<i4'), ('galid', '<i4'), ('hhost', '<i4', (5,)),
                    ('ang', '<f4', (3,)), ('m', '<f4'), ('mvir', '<f4'),
                    ('r', '<f4', (4,)), ('energy', '<f8', (4,)),
                    ('rvir', '<f4'), ('p', '<f4', (3,)),
                    ('v', '<f4', (3,)), ('sp', '<f4')]
            
            dtype_profile = [('rho_0', '<f4'), ('r_c', '<f4')]
            
            if tothal is None:
                tothal = self.halnum + self.subnum
            self.data = np.recarray(tothal, dtype=dtype_gal)
            self.profile = np.recarray(tothal, dtype=dtype_profile)
            for i in range(5):
                nph = read_fortran(f, np.dtype('i4'), 1)
                self.data['np'][i] = nph
                pids = read_fortran(f, np.dtype('i4'), nph) # id list. 
                self.data['galid'][i] = read_fortran(f, np.dtype('i4'), 1)
                read_fortran(f, np.dtype('i4'), 1) #timestep
                self.data['hhost'][i][0:5] = read_fortran(f, np.dtype('i4'), 5)
                self.data['m'][i] = read_fortran(f, np.dtype('f4'), 1)
                self.data['p'][i][0:3] = read_fortran(f, np.dtype('f4'), 3)
                self.data['v'][i][0:3] = read_fortran(f, np.dtype('f4'), 3)
                self.data['ang'][i][0:3] = read_fortran(f, np.dtype('f4'), 3)
                self.data['r'][i][0:4] = read_fortran(f, np.dtype('f4'), 4)
                self.data['energy'][i][0:4] = read_fortran(f, np.dtype('f4'), 3)#energies
                self.data['sp'][i] = read_fortran(f, np.dtype('f4'), 1)
                sig, sigb, mb = read_fortran(f, np.dtype('f4'), 3)
                virials = read_fortran(f, np.dtype('f4'), 4)
                self.data['rvir'][i], self.data['mvir'][i] = virials[0:2]
                self.profile['rho_0'][i], self.profile['r_c'] = read_fortran(f, np.dtype('f4'), 2)
                nbin = read_fortran(f, np.dtype('i4'), 1)
                rr = read_fortran(f, np.dtype('f4'), nbin)
                rho = read_fortran(f, np.dtype('f4'), nbin) # wasted..
                
               
                
                
                
                
                
                
                