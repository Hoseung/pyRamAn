# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:32:55 2015

@author: hoseung
"""
import numpy as np
from load.utils import read_header, read_fortran, skip_fortran
from load.sim import Simbase

class Dummy():
    def __init__(self):
        pass


class Hydro(Simbase):
    """
    
    Notes
    -----
    Hydro, part, amr share some methods and attributes such as set_info, set_ranges, out_dir, _fnbase, and so on. 
    I could make a base class for raw data with the common methods and attributes,
    and inherite the base class to define Hydro, Part, and Amr class.
    """

    def __init__(self, nout=None, info=None, amr=None, region=None, ranges=None, load=False):
        """
        Parameters
        ----------
        region : dict-like
            d
        ranges : array-like (3 by 2)
            region preceeds(?) ranges.
        """
        if info is None:
            assert nout is not None, "either info or onut is required"
            from load.info import Info
            print("[Hydro.__init__] Loading info")
            info = Info(nout=nout)
        self.info = info
        self.nout = info.nout

        snout = str(info.nout).zfill(5)
        # file name
        self.out_dir = 'snapshots/'
        self._fnbase = info.base + '/' + self.out_dir + 'output_' + snout + '/hydro_' + snout + '.out'
        self._get_basic_info()
        self.set_info(info)
        if region is not None:
            ranges = region['ranges']
        if ranges is not None:
            self.set_ranges(ranges=ranges)
        elif info.ranges is not None:
            self.set_ranges(ranges=info.ranges)

        try:
            self.amr = amr
        except NameError:
            print("Loading amr first! \n")

        if load:
            self.amr2cell()

    def _get_basic_info(self):
        f = open(self._fnbase + '00001', "rb")
        self.header = Dummy()
        self._read_hydro_header(f)

    def set_info(self, info):
        self.info = info

    def _read_hydro_header(self, f, verbose=False):
        # Global
        h1 = read_header(f,
                         np.dtype(
                             [('ncpu', 'i'),
                              ('nvarh', 'i4'),
                              ('ndim', 'i4'),
                              ('nlevelmax', 'i4'),
                              ('nboundary', 'i4'),
                              ('gamma', 'f8')]))

        if h1['nboundary'] == 0 and verbose:
            print(' Periodic boundary condition')

        # if assign
        self.header.ncpu = h1['ncpu']
        self.header.ndim = h1['ndim']
        self.header.nlevelmax = h1['nlevelmax']
        self.header.nboundary = h1['nboundary']
        self.header.gamma = h1['gamma']
        self.header.nvarh = h1['nvarh']

    def amr2cell(self, lmax=None, icpu=0, cpu=True,
                 verbose=False, return_meta=False,
                 ranges=None):
        """
        Loads AMR and HYDRO and output hydro data into particle-like format(cell).

        Parameters
        ----------
        cpu : bool, optional
            If True, cpu number of each cell is stored.
        icpu : int, array-like, optional
            list of cpus to load, has no effect... 
        lmax : int, optional
            Limit the maximum level of hydro data returned.
        return_meta : bool, optional
            If True, returns meta data instead. (Why would I want that??)
        verbose : bool, optional
            
        """
        from load import a2c
        nvarh = self.header.nvarh
        nlevelmax = self.header.nlevelmax
        if lmax is None:
            lmax = nlevelmax

        if verbose: print(' >>> working resolution (lmax) =', lmax)

        if ranges is not None: self.set_ranges(ranges=ranges)
        # Set ranges
        xmi, xma = self.ranges[0]
        ymi, yma = self.ranges[1]
        zmi, zma = self.ranges[2]

        work_dir = self.info.base + '/' + self.out_dir + 'output_' + str(self.info.nout).zfill(5)
        if verbose:
            print("[hydro.amr2cell] Ranges", xmi, xma, ymi, yma, zmi,zma)
            print("[hydro.amr2cell] cpus", self.cpus)
            
        out = a2c.a2c_count(work_dir, xmi, xma, ymi, yma, zmi, zma, lmax, self.cpus)

        if return_meta:
            return (out[0], work_dir, xmi, xma, ymi, yma, zmi, zma, lmax)
        else:
            cell = a2c.a2c_load(work_dir, xmi, xma, ymi, yma, zmi, zma,\
                                lmax, out[0], nvarh, self.cpus)
            dtype_cell = [('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('dx', '<f8')] 
            if cpu:
                dtype_cell.append(('cpu', '<f8'))
            for i in range(nvarh):
                dtype_cell.append( ('var' + str(i), '<f8'))

        self.cell = np.zeros(len(cell[1]), dtype=dtype_cell)
        self.cell['x'] = cell[0][:,0]
        self.cell['y'] = cell[0][:,1]
        self.cell['z'] = cell[0][:,2]
        self.cell['dx'] = cell[1]
        for i in range(nvarh):
            self.cell['var' + str(i)] = cell[2][:,i]
        if cpu:
            self.cell['cpu'] = cell[3]
#        self.cell = np.rec.fromarrays([cell[0][:,0], cell[0][:,1], cell[0][:,2], cell[1],
#                                       cell[2][:,0], cell[2][:,1], cell[2][:,2],
#                                       cell[2][:,3], cell[2][:,4], cell[2][:,5]],
#                                       dtype = dtype_cell)

