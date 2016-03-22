# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:32:55 2015

@author: hoseung
"""
import numpy as np
from load.utils import read_header, read_fortran, skip_fortran
import load

class Hydro:
    class HydroHeader():
        def __init__(self):
            pass

    def __init__(self, info, amr=None):

        snout = str(info.nout).zfill(5)
        # file name
        self._fnbase = info.base + '/snapshots/output_' + snout + '/hydro_' + snout + '.out'
        self._get_basic_info()
        self.set_info(info)
        try:
            self.amr = amr
        except NameError:
            print("Load amr first! \n")

    def _get_basic_info(self):
        f = open(self._fnbase + '00001', "rb")
        self.header = self.HydroHeader()
        self._read_hydro_header(f)

    def set_info(self, info):
        self.info = info

    def _read_hydro_header(self, f, verbose=False):
        # Global
        h1 = read_header(f,
                         np.dtype(
                             [('ncpu', 'i4'),
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

    def load(self, lmax=None, icpu=0):
        print("load")
        
    def amr2cell(self, lmax=None, icpu=0, verbose=False, return_meta=False):
        from load import a2c
        # If there's no source code, compile it with scipy. 
#        import numpy as np
#        cpus = self.info.cpus
#        ndim = self.info.ndim
#        print("CPU list:", cpus)
        nvarh = self.header.nvarh
#        ncpu = self.header.ncpu
#        nboundary = self.header.nboundary
        nlevelmax = self.header.nlevelmax
        if lmax is None:
            lmax = nlevelmax

        print(' >>> working resolution (lmax) =', lmax)
        #print("a2ca2c")

        # Set ranges
        xmi = self.info.ranges[0][0]
        xma = self.info.ranges[0][1]
        ymi = self.info.ranges[1][0]
        yma = self.info.ranges[1][1]
        zmi = self.info.ranges[2][0]
        zma = self.info.ranges[2][1]
        work_dir = self.info.base + '/snapshots/output_' + str(self.info.nout).zfill(5)
        
        out = a2c.a2c_count(work_dir, xmi, xma, ymi, yma, zmi, zma, lmax)
        if return_meta:
            return (out[0], work_dir, xmi, xma, ymi, yma, zmi, zma, lmax)

        else:
            cell = a2c.a2c_load(work_dir, xmi, xma, ymi, yma, zmi, zma, lmax
                                 ,out[0], nvarh)
            dtype_cell = [('x', '<f8'), ('y', '<f8'), ('z', '<f8'),('dx', '<f8'),
                      ('var0', '<f8'), ('var1', '<f8'), ('var2', '<f8'),
                        ('var3', '<f8'), ('var4', '<f8'), ('var5', '<f8')]

        self.cell = np.rec.fromarrays([cell[0][:,0], cell[0][:,1], cell[0][:,2], cell[1],
                                       cell[2][:,0], cell[2][:,1], cell[2][:,2],
                                       cell[2][:,3], cell[2][:,4], cell[2][:,5]],
                                       dtype = dtype_cell)


    def amr2cell_old(self, lmax=None, icpu=0, verbose=False):
        """
        Reads AMR and HYDRO and output hydro variables in particle-like data format.
        Works ok.

        Initial version: 2015. 04. 16
        """
        import numpy as np
        # if both amr and hydro does not exist, create.
        #################
        cpus = self.info.cpus
        ndim = self.info.ndim
        print("CPU list:", cpus)
        twotondim = 2**ndim
        nvarh = self.header.nvarh
        ncpu = self.header.ncpu
        nboundary = self.header.nboundary
        nlevelmax = self.header.nlevelmax
        if lmax is None:
            lmax = nlevelmax

        print(' >>> working resolution (lmax) =', lmax)

        # Set ranges
        xmi = self.info.ranges[0][0]
        xma = self.info.ranges[0][1]
        ymi = self.info.ranges[1][0]
        yma = self.info.ranges[1][1]
        zmi = self.info.ranges[2][0]
        zma = self.info.ranges[2][1]

        # Size of arrays
        if icpu == 0:
            listmax = ncpu
        elif icpu > 0:
            listmax = ncpu + nboundary
        ngridarr = np.zeros((nlevelmax, listmax), dtype=np.int32)

        nx = self.amr.header.ng[0]
        ny = self.amr.header.ng[1]
        nz = self.amr.header.ng[2]
        xbound = [int(nx/2), int(ny/2), int(nz/2)]

        # Initialize pointer array

        nlinemax = 10000000
        xt = np.zeros(nlinemax, dtype=np.float64)
        yt = np.zeros(nlinemax, dtype=np.float64)
        zt = np.zeros(nlinemax, dtype=np.float64)
        dxt = np.zeros(nlinemax, dtype=np.float64)
        vart = np.zeros((nvarh, nlinemax), dtype=np.float64)

        nline = 0
        for icpu in cpus:
            print("reading", icpu)
            # Open AMR file and read header
            famr = open(self.amr._fnbase + str(icpu).zfill(5), "rb")

            header_icpu = load.amr.AmrHeader()
            header_icpu._read_amr_header(famr)

            # This is important ??
            ngridarr[0: nlevelmax][0: ncpu + nboundary] = header_icpu.numbl
#            if nboundary > 0:
#          ngridarr[0: nlevelmax][ncpu: ncpu + nboundary] = header_icpu.numbb

            # because hydro file is written in fortran order
            # looping over
            # for ilevel
            #   for icpu
            # Open Hydro file and skip header

            fhydro = open(self._fnbase + str(icpu).zfill(5), "rb")
            self._read_hydro_header(fhydro)

            xc = np.zeros((3, 8), dtype=np.float64)

            for ilevel in range(1, lmax + 1):
                # geometry
                dx = 0.5**ilevel  # grid size of the current level.
                dx2 = dx/2  # half of the grid size of the current level.
                nx_full = 2**ilevel  # maximum possible number of grid cells at the given level.
                ny_full = 2**ilevel
                nz_full = 2**ilevel
                # division results in floats!
                for ioct in range(twotondim):  #  1 to 8
                    iz = int((ioct) / 4)
                    iy = int((ioct - 4 * iz) / 2)
                    ix = int((ioct - 2 * iy - 4 * iz))
                    xc[0][ioct] = (ix - 0.5) * dx  # grid center
                    xc[1][ioct] = (iy - 0.5) * dx
                    xc[2][ioct] = (iz - 0.5) * dx

                # Allocate work arrays
                ngrida = ngridarr[ilevel-1][icpu-1]

                if ngrida > 0:
                    # note that index stars from 0 to reduce memory use
                    xg = np.zeros((ndim, ngrida), dtype=np.float64)
                    son = np.zeros((twotondim, ngrida), dtype=np.int32)
                    var = np.zeros((nvarh, twotondim, ngrida,),
                                   dtype=np.float64)
                    x = np.zeros((ndim, ngrida), dtype=np.float64)
#                    rho = np.zeros(ngrida, dtype=np.float64)

# Main loop, loop over doamin
                for jcpu in range(1, nboundary + ncpu + 1):
                    ng = ngridarr[ilevel-1][jcpu-1]
                    if ng > 0: # If ng == 0, go straight to Read Hydro data
                        if (verbose):
                            print("Level {:2d} has {:6d} grids in proc {:4d}".format(ilevel, ng, jcpu))
# SKIP
                        read_fortran(famr, np.int32, ng)  # grid index
                        read_fortran(famr, np.int32, ng)  # next index
                        read_fortran(famr, np.int32, ng)  # prev index
# Read gird center
                        for idim in range(ndim):
                            if jcpu == icpu:
                                xtemp = read_fortran(famr, np.dtype('f8'), ng)
                                xg[idim][:] = xtemp
                            else:
                                read_fortran(famr, np.dtype('f8'), ng)
# skip father index
                        read_fortran(famr, np.int32, ng)
# skip neighbour index ( 6 possible neighbours )
                        for idim in range(2*ndim):
                            read_fortran(famr, np.int32, ng)
# Read son index ( 8 possible children, oct-tree )
                        for idim in range(2**ndim):
                            if jcpu == icpu:
                                son[idim][:] = read_fortran(famr, np.int32, ng)
                            else:
                                read_fortran(famr, np.int32, ng)
# skip cpu map
                        for idim in range(2**ndim):
                            read_fortran(famr, np.dtype('i4'), ng)
# skip ref map
                        for idim in range(2**ndim):
                            read_fortran(famr, np.dtype('i4'), ng)
# Read Hydro data
# Skip Hydro header always.
                    skip_fortran(fhydro)
                    skip_fortran(fhydro)

                    if ng > 0:
                        for ioct in range(twotondim):
                            for ivar in range(nvarh):
                                if jcpu == icpu:
                                    var[ivar][ioct][:] = read_fortran(fhydro, np.dtype('f8'), ng)
                                else:
                                    read_fortran(fhydro, np.dtype('f8'), ng)
# Compute map
                if ngrida > 0:
                    for ioct in range(twotondim):
                        for idim in range(3):
                            x[idim][:] = xg[idim][:] + xc[idim][ioct] - xbound[idim]

                        ref = np.zeros(ngrida, dtype=np.int16)
                        if ilevel < lmax:
                            val = np.where(son[ioct][:] > 0)[0]
                            if len(val) > 0:
                                ref[val] = 1
                        ok_cell = np.where((ref == 0) &
                                           (x[0, :] + dx2 >= xmi) &
                                           (x[1, :] + dx2 >= ymi) &
                                           (x[2, :] + dx2 >= zmi) &
                                           (x[0, :] - dx2 <= xma) &
                                           (x[1, :] - dx2 <= yma) &
                                           (x[2, :] - dx2 <= zma))[0]

#################  It is extremely slow!!!   Use np.where instead.
#                            ok_cell = np.all([ref == 0,
#                                              x[0, :] + dx2 >= xmi, x[0, :] - dx2 <= xma,
#                                              x[1, :] + dx2 >= ymi, x[1, :] - dx2 <= yma,
#                                              x[2, :] + dx2 >= zmi, x[2, :] - dx2 <= zma],
#                                              axis=0)

                        Nok_cell = len(ok_cell)
                        if Nok_cell > 0:
                            nbegin = nline
                            nend = nline + Nok_cell

                            # If the existing arrays are not large enough.
                            if nend > nlinemax:
                                print(nend, " > ", nlinemax)
                                xt = np.append(xt, np.zeros(nlinemax, dtype=np.float64))
                                yt = np.append(yt, np.zeros(nlinemax, dtype=np.float64))
                                zt = np.append(zt, np.zeros(nlinemax, dtype=np.float64))
                                dxt = np.append(dxt, np.zeros(nlinemax, dtype=np.float64))
                                vart = np.append(vart, np.zeros((nvarh,nlinemax), dtype=np.float64), axis=1)
                                nlinemax += nlinemax
# save variables
                            xt[nbegin:nend] = x[0, ok_cell]
                            yt[nbegin:nend] = x[1, ok_cell]
                            zt[nbegin:nend] = x[2, ok_cell]
                            dxt[nbegin:nend] = dx
                            vart[:, nbegin:nend] = var[:, ioct, ok_cell]
                            nline = nend

            # loop over ilevel
            famr.close()
            fhydro.close()

        dtype_cell = [('x', '<f8'), ('y', '<f8'), ('z', '<f8'),('dx', '<f8'),
                      ('var0', '<f8'), ('var1', '<f8'), ('var2', '<f8'),
                        ('var3', '<f8'), ('var4', '<f8'), ('var5', '<f8')]

        self.cell = np.rec.fromarrays([xt[:nend], yt[:nend], zt[:nend], dxt[:nend],
                                       vart[0][:nend], vart[1][:nend], vart[2][:nend],
                                       vart[3][:nend], vart[4][:nend], vart[5][:nend]],
                                       dtype = dtype_cell)

    def passpass(self):
        print("   ")
