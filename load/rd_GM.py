# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 04:50:37 2016

@author: hoseung
"""
from load.utils import read_header, read_fortran
import numpy as np


class Gal():
    def __init__(self):
        self.data = None
        self.header = None


def rd_gal(nout, idgal, wdir="./", metal=True, nchem=0,long=True):
    # Header structure
    dtype_header = np.dtype([('my_number', 'i4'),
                             ('level', 'i4'),
                             ('mgal', 'f8'),
                             ('xg', 'f8', (3,)),
                             ('vg', 'f8', (3,)),
                             ('lg', 'f8', (3,)),
                             ('npart', 'i4')])
    
    # variable type
    dtype_data=[('pos', '<f8', (3,)),
           ('vel', '<f8', (3,)),
           ('id', '<i4'),
           ('m', '<f8'),
           ('time', '<f8')]
    if metal:
        dtype_data.append(('metal', '<f8'))
        if nchem > 0:
            dtype_data.append(('cp', '<f8', (nchem,)))
    
    idgal = str(idgal).zfill(7)
    dir_nout = "GAL_" + str(nout).zfill(5)      
    with open(wdir + 'GalaxyMaker/' +  dir_nout + '/gal_stars_' + idgal, "rb") as f:
        header = read_header(f, dtype=dtype_header)
        header['mgal'] *= 1e11 # mass fof galaxy
        
        # data array
        data = np.recarray(header["npart"], dtype=dtype_data)
        data['pos'][:,0] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['pos'][:,1] = read_fortran(f, np.dtype('f8'), header["npart"])
        data['pos'][:,2] = read_fortran(f, np.dtype('f8'), header["npart"])
    
        data['vel'][:,0] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['vel'][:,1] = read_fortran(f, np.dtype('f8'), header["npart"])
        data['vel'][:,2] = read_fortran(f, np.dtype('f8'), header["npart"])
        
        
        data['m'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['id'] = read_fortran(f, np.dtype('i4'), header["npart"])
        data['time'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        if 'metal' in data.dtype.names:
            data['metal'] = read_fortran(f, np.dtype('f8'), header["npart"])
            if nchem > 0:
                for i in range(nchem):
                    data['cp'][:,i] = read_fortran(f, np.dtype('f8'), header["npart"])
        
    gal = Gal()
    gal.data = data
    gal.header = header
    return gal
