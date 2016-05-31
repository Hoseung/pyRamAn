# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 04:50:37 2016

@author: hoseung
"""
from load.utils import read_header, read_fortran
import numpy as np


class Gal():
    def __init__(self):
        """
        
        """
        self.star = None # data -> star
        self.header = None


def rd_gal(nout, idgal, wdir="./", metal=True, nchem=0,long=True):
    """

    header xg in Mpc (physical, centered at 0.5, 0.5, 0.5 of the simualtion volume)
    """
    # Header structure
    dtype_header = np.dtype([('my_number', 'i4'),
                             ('level', 'i4'),
                             ('mgal', 'f8'),
                             ('xg', 'f8', (3,)),
                             ('vg', 'f8', (3,)),
                             ('lg', 'f8', (3,)),
                             ('npart', 'i4')])
    
    # variable type
    dtype_data=[('x', '<f8'),
                ('y', '<f8'),
                ('z', '<f8'),
                ('vx', '<f8'),
                ('vy', '<f8'),
                ('vz', '<f8'),
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
        data['x'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['y'] = read_fortran(f, np.dtype('f8'), header["npart"])
        data['z'] = read_fortran(f, np.dtype('f8'), header["npart"])
    
        data['vx'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['vy'] = read_fortran(f, np.dtype('f8'), header["npart"])
        data['vz'] = read_fortran(f, np.dtype('f8'), header["npart"])
        
        
        data['m'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['id'] = read_fortran(f, np.dtype('i4'), header["npart"])
        data['time'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        if 'metal' in data.dtype.names:
            data['metal'] = read_fortran(f, np.dtype('f8'), header["npart"])
            if nchem > 0:
                for i in range(nchem):
                    data['cp'][:,i] = read_fortran(f, np.dtype('f8'), header["npart"])
        
    gal = Gal()
    gal.star = data
    gal.header = header
    return gal


def rd_cell(nout, idgal, wdir="./", metal=True, nchem=0):
    from utils import io
    dir_nout = "CELL_" + str(nout).zfill(5)   
    with open(wdir + 'GalaxyMaker/' +  dir_nout + '/gal_cells_'+ 
    str(idgal).zfill(7), "rb") as f:
        nout0 = io.read_fortran(f, dtype=np.int32, check=False)[0]
        assert nout == nout0, "given nout ({}) and loaded nout ({}) do not match".format(nout, nout0)
        gid = io.read_fortran(f, dtype=np.int32, check=False)[0] 
        assert idgal == gid, "given idgal ({}) and loaded idgal ({}) do not match".format(idgal, gid)
        
        ncell = io.read_fortran(f, dtype=np.int32, check=False)[0]
        cell = np.zeros(ncell, dtype=[('x', '<f8'),('y', '<f8'),('z', '<f8'),
                                     ('dx', '<f8'),('var0', '<f8'),('var1', '<f8'),
                                     ('var2', '<f8'),('var3', '<f8'),('var4', '<f8'),
                                     ('var5', '<f8')])
        cell['x']  = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['y']  = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['z']  = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['dx'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var0'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var1'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var2'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var3'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var4'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var5'] = io.read_fortran(f, dtype=np.float64, n=ncell)
    return cell
