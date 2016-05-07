# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 14:02:51 2016

@author: hoseung
"""

def load_meta(f):
    import numpy as np
    from load.utils import read_fortran
    read_fortran(f, np.dtype('i4'), 1)
    read_fortran(f, np.dtype('f4'), 1)
    read_fortran(f, np.dtype('f4'), 1)
    read_fortran(f, np.dtype('f4'), 1)
    read_fortran(f, np.dtype('f4'), 1)
    halnum, subnum = read_fortran(f, np.dtype('i4'), 2)
    return halnum, subnum


def load_data(f, tothal):
    import numpy as np
    from load.utils import read_fortran
    dtype_halo = [('np', '<i4'), ('hnu', '<i4'), ('hhost', '<i4', (5,)),
          ('ang', '<f8', (3,)), ('m', '<f8'), ('mvir', '<f8'),
            ('r', '<f8'), ('rvir', '<f8'), ('p', '<f8', (3,)),
            ('v', '<f8', (3,)), ('sp', '<f8')]

    data = np.recarray(tothal, dtype=dtype_halo)
    for i in range(tothal):
        nph = read_fortran(f, np.dtype('i4'), 1)
        data['np'][i] = nph
        read_fortran(f, np.dtype('i4'), nph) # id list. 
        hnu = read_fortran(f, np.dtype('i4'), 1)
        data['hnu'][i] = hnu
        
        read_fortran(f, np.dtype('i4'), 1) #timestep
        data['hhost'][i][0:5] = read_fortran(f, np.dtype('i4'), 5)
        data['m'][i] = read_fortran(f, np.dtype('f4'), 1)
        data['p'][i][0:3] = read_fortran(f, np.dtype('f4'), 3)
        data['v'][i][0:3] = read_fortran(f, np.dtype('f4'), 3)
        data['ang'][i][0:3] = read_fortran(f, np.dtype('f4'), 3)
        data['r'][i] = read_fortran(f, np.dtype('f4'), 4)[0]
        read_fortran(f, np.dtype('f4'), 3)#energies
        data['sp'][i] = read_fortran(f, np.dtype('f4'), 1)
        data['rvir'][i], data['mvir'][i] = read_fortran(f, np.dtype('f4'), 4)[0:2]
        read_fortran(f, np.dtype('f4'), 2)
        
        
if __name__ == "__main__":
    import time
    f = open("~/Work/data/01605/halo/DM/tree_bricks187", "rb")
    numh, nums = load_meta(f)
    t = time.time()
    load_data(f, numh + nums)
    print("Took {} seconds".format(time.time() - t))
    