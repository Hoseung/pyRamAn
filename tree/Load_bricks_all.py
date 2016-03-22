# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 09:22:32 2015

@author: hoseung
"""

def load_brick_all(f):
    import numpy as np
    from load.utils import read_fortran, skip_fortran
    
    is_gal=True
   
    
    nbodies = read_fortran(f, np.dtype('i4'), 1)
    massp = read_fortran(f, np.dtype('f4'), 1)
    aexp = read_fortran(f, np.dtype('f4'), 1)
    omegat = read_fortran(f, np.dtype('f4'), 1)
    age = read_fortran(f, np.dtype('f4'), 1)
    halnum, subnum = read_fortran(f, np.dtype('i4'), 2)        
    
    dtype_halo = [('np', '<i4'), ('id', '<i4'), ('level', '<i4'),
                  ('host', '<i4'), ('sub', '<i4'), ('nsub', '<i4'),
                ('nextsub', '<i4'),
                ('m', '<f4'), ('mvir', '<f4'),
                ('r', '<f4'), ('rvir', '<f4'), 
                ('x', '<f4'), ('y', '<f4'), ('z', '<f4'),
                ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4'),
                ('ax', '<f4'), ('ay', '<f4'), ('az', '<f4'),
                ('sp', '<f4')]
    
    tothal = halnum + subnum
    
    data = np.recarray(tothal, dtype=dtype_halo)
    idlists=[[]]*tothal
    for i in range(tothal):
        nph = read_fortran(f, np.dtype('i4'), 1)
        data['np'][i] = nph
        idlists[i] = read_fortran(f, np.dtype('i4'), nph) # id list. 
        data['id'][i] = read_fortran(f, np.dtype('i4'), 1)
        read_fortran(f, np.dtype('i4'), 1) #timestep
        data['level'][i], data['host'][i], \
        data['sub'][i], data['nsub'][i], \
        data['nextsub'][i] = read_fortran(f, np.dtype('i4'), 5)
        data['m'][i] = read_fortran(f, np.dtype('f4'), 1)
        data['x'][i], data['y'][i], data['z'][i] \
        = read_fortran(f, np.dtype('f4'), 3)
        data['vx'][i], data['vy'][i], data['vz'][i] \
        = read_fortran(f, np.dtype('f4'), 3)
        data['ax'][i], data['ay'][i], data['az'][i] \
        = read_fortran(f, np.dtype('f4'), 3)                
        data['r'][i] = read_fortran(f, np.dtype('f4'), 4)[0]
        read_fortran(f, np.dtype('f4'), 3)#energies
        data['sp'][i] = read_fortran(f, np.dtype('f4'), 1)
        if is_gal:
            #self.data['sig'][i], self.data['sigb'][i], self.data['mb'][i] \
            skip_fortran(f)
            #= read_fortran(f, np.dtype('f4'), 3)
        data['rvir'][i], data['mvir'][i] = read_fortran(f, np.dtype('f4'), 4)[0:2]
        read_fortran(f, np.dtype('f4'), 2)
        if is_gal:
            skip_fortran(f)
            skip_fortran(f)
            skip_fortran(f)
        
    f.close()
return data, idlists

#%%


f = open('/home/hoseung/Work/data/05427/GalaxyMaker/gal/tree_bricks186', 'rb')
data0, idlists0 = load_brick_all(f)
f = open('/home/hoseung/Work/data/05427/GalaxyMaker/gal/tree_bricks187', 'rb')
data1, idlists1 = load_brick_all(f)

#%%