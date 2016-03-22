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
        self.info = info
        try:
            self.amr = amr
        except NameError:
            print("Load amr first! \n")

    def _get_basic_info(self):
        f = open(self._fnbase + '00001', "rb")
        self.header = self.HydroHeader()
        self._read_hydro_header(f)

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

    def amr2cell(self, lmax=None, icpu=0, verbose=False):
        # if both amr and hydro does not exist, create.
        #################
        cpus = self.info.cpus
        ndim = self.info.ndim
        print("CPU list:", cpus)
        twotondim = 2**ndim
#        ngridtot = 0

        nvarh = self.header.nvarh
        ncpu = self.header.ncpu
        nboundary = self.header.nboundary
        nlevelmax = self.header.nlevelmax
        if lmax is None:
            lmax = nlevelmax

        print(' >>> working resolution (lmax) =', lmax)

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
        xbound = [nx/2, ny/2, nz/2]
        print("xbound:",xbound)

        nlinemax = 20000000
        nline = 0
        for icpu in cpus:
            # Open AMR file and read header
            famr = open(self.amr._fnbase + str(icpu).zfill(5), "rb")

            header_icpu = load.amr.AmrHeader()
            header_icpu._read_amr_header(famr)

            # This is important

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

            # Initialize pointer array
            xg = np.zeros((ndim, nlinemax), dtype=np.float64)
            xt = np.zeros(nlinemax, dtype=np.float64)
            yt = np.zeros(nlinemax, dtype=np.float64)
            zt = np.zeros(nlinemax, dtype=np.float64)
            dxt = np.zeros(nlinemax, dtype=np.float64)
            vart = np.zeros((nvarh, nlinemax), dtype=np.float64)

            xc = np.zeros((3, 8), dtype=np.float64)

            for ilevel in range(1, lmax + 1):
                # geometry
                dx = 0.5**ilevel  # grid size of the current level.
                dx2 = dx/2  # half of the grid size of the current level.
#                print('level',ilevel,'dx2',dx2)
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

    #           grid[ilevel].ngrid=ngrida
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
                    if ng > 0:
                        if (verbose):
                            print("Level %2d has %6d grids in proc %4d"
                                  % (ilevel, ng, jcpu))
# SKIP
                        read_fortran(famr, np.int32, ng)  # grid index
                        read_fortran(famr, np.int32, ng)  # next index
                        read_fortran(famr, np.int32, ng)  # prev index
# Read gird center
                        for idim in range(ndim):
                            if jcpu == icpu:
                                xg[idim][:] = read_fortran(famr, np.dtype('f8'), ng)
                            else:
                                read_fortran(famr, np.dtype('f8'), ng)
# skip father index
                        read_fortran(famr, np.int32, ng)
# skip neighbour index ( 6 possible neighbours)
                        for idim in range(2*ndim):
                            read_fortran(famr, np.int32, ng)
# skip son index ( 8 possible children, oct-tree)
                        for idim in range(2**ndim):
                            temp = read_fortran(famr, np.int32, ng)
                            if jcpu == icpu:
                                son[idim][:] = temp
# skip cpu map
                        for idim in range(2**ndim):
                            read_fortran(famr, np.dtype('i4'), ng)
# skip ref map
                        for idim in range(2**ndim):
                            read_fortran(famr, np.dtype('i4'), ng)

# Read Hydro data #
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
#                            x[0][:] = xg[:][0] + xc[0][ioct]-xbound[0]
#                            x[1][:] = xg[:][1] + xc[1][ioct]-xbound[1]
#                            x[2][:] = xg[:][2] + xc[2][ioct]-xbound[2]

                            ref = np.zeros(ngrida, dtype=np.int16)
                            if ilevel < lmax:
                                val = np.where(son[ioct][:] > 0)[0]
                                if len(val) > 0:
#                                    print(ngrida,'aaa',len(val))
                                    ref[val] = 1

#                            ok_cell = np.where((ref == 0) &
#                                               (x[0, :] + dx2 >= xmi))# &
#                                               (x[1, :] + dx2 >= ymi) &
#                                               (x[2, :] + dx2 >= zmi) &
#                                               (x[0, :] - dx2 <= xma) &
#                                               (x[1, :] - dx2 <= yma) &
#                                               (x[2, :] - dx2 <= zma))

                            ok_cell = np.all([ref == 0,
                                              x[0, :] + dx2 >= xmi, x[0, :] - dx2 <= xma,
                                              x[1, :] + dx2 >= ymi, x[1, :] - dx2 <= yma,
                                              x[2, :] + dx2 >= zmi, x[2, :] - dx2 <= zma],
                                              axis=0)

                            Nok_cell = sum(ok_cell)
#                            Nok_cell = len(ok_cell[0])

                            if Nok_cell > 0:
                                print('ok CELL', Nok_cell)
                                nbegin = nline
                                nend = nline + Nok_cell - 1

                                # If the existing arrays are not large enough.
                                # Must be better ways.
                                if nend > nlinemax:
                                    xt = np.append(xt, np.zeros(nlinemax, dtype=np.float64))
                                    yt = np.append(yt, np.zeros(nlinemax, dtype=np.float64))
                                    zt = np.append(zt, np.zeros(nlinemax, dtype=np.float64))
                                    dxt = np.append(dxt, np.zeros(nlinemax, dtype=np.float64))
                                    vart = np.append(vart, np.zeros((nvarh,nlinemax), dtype=np.float64))

                                xt[nbegin:nend] = x[0, ok_cell]
                                yt[nbegin:nend] = x[1, ok_cell]
                                zt[nbegin:nend] = x[2, ok_cell]
                                dxt[nbegin:nend] = dx[ok_cell]
                                vart[:, nbegin:nend] = var[:, ioct, ok_cell]

                                nline = nend+1
                                print(nline)

            # loop over ilevel
            famr.close()
            fhydro.close()

        print(xt,yt,vart[0])
        dtype_cell = [('x', '<f8'), ('y', '<f8'), ('z', '<f8'),
                      ('dx', '<f8'), ('var1', '<f8'),
                      ('var2', '<f8'), ('var3', '<f8'), ('var4', '<f8'),
                      ('var5', '<f8'), ('var6', '<f8')]
        self.cell = np.rec.fromarrays([xt, yt, zt, dxt,
                                       vart[0], vart[1], vart[2],
                                       vart[3], vart[4], vart[5]],
                                       dtype = dtype_cell)
                                       
                                       
#%%
from load import amr
from load import info
from load import part
class Sim():
    """
    Defines the 'host class' of part, amr, hydro, info.

    The most global information are stored directly to this class:
    ndim, ncpu, base.
    Later it will also include .nml information.
    (Roman's git version generates such output in text files)

    Currently this class deals with single snapshot.
    But I hope to expand it for multiple snapshots.
    """
    def __init__(self, nout, base='./', ranges=[[0.0,1.0],[0.0,1.0],[0.0,1.0]]):
        self.nout = nout
        self.set_base(base)
        self.add_info()
        self.ranges = []
        self.set_ranges(ranges)
        self.add_amr()
        self.set_cpus(self._hilbert_cpulist(self.info, self.ranges))
        print(' ')

    def add_info(self):        
#        from load import info
        self.info = info.Info(self.nout, self.base)

    def add_part(self, ptypes):
        print("Types of particles you want to load are: ", ptypes)
        self.part = part.Part(self.info, ptypes)
        print("A particle instance is created")
        print("Use part.load() to load these particle")

    def add_amr(self):
        self.amr = amr.Amr(self.info)
        print("An AMR instance is created\n")

    def add_hydro(self):
        self.hydro = Hydro(self.info, self.amr) ###
        print("An Hydro instance is created\n")

    def set_base(self, base):
        """
            Set Working directory.
        """
        self.base = base
        self.show_base()

    def show_base(self):
        print("setting the base(working) directory to :", self.base)

    def set_cpus(self, cpus):
        self.cpus = cpus
        try:
            print("Updating info.cpus")
            self.info._set_cpus(self.cpus)
        except AttributeError:
            print("No info._set_cpus attribute??")
        self.show_cpus()

    def show_cpus(self):
        print(" ncpus : %s \n" % self.cpus)

    def set_ranges(self, ranges=[[0, 1], [0, 1], [0, 1]]):
        nr = np.asarray(ranges)
        if not(nr.shape[0] == 3 and nr.shape[1] == 2):
            print('Shape of ranges is wrong:', nr.shape)
            print('example : [[0.1,0.3],[0.2,0.4],[0.6,0.8]]')
        else:
            self.ranges = ranges
            try:
                self.info._set_ranges(self.ranges)
            except AttributeError:
                print("There is no info._set_ranges attribute")
            self.show_ranges()
            self.set_cpus(self._hilbert_cpulist(self.info, self.ranges))

    def show_ranges(self):
        print("Ranges = %s\n" % self.ranges)

    def zoomin_region(self, *args, **kwargs):  # If part is not loaded yet, load particles
        """
        Not only part, but also hydro or amr can be used to find the zoomin region!
        """
        self.part.get_zoomin( *args, **kwargs)

    def _hilbert_cpulist(self, info, ranges):
        '''
        After determining cpu numbers, read the cpu files and cut off data
        that are outside ranges.
        -> cpu files contain data points within a given ragnes
        BUT! they also contain other data points outside the ragnes.
        '''
        if not(hasattr(self, 'amr')):
            self.add_amr()
        nlevelmax = self.amr.header.nlevelmax
        nboundary = self.amr.header.nboundary
        ndim = self.info.ndim
        ncpu = self.info.ncpu_tot  #ncpu_tot
        # header of part output is needed (and is prepared.)

        # header of amr output is needed.
        two2ndim = 2**ndim
        nx = self.amr.header.ng[0]
        ny = self.amr.header.ng[1]
        nz = self.amr.header.ng[2]

        xbound = [nx/2., ny/2., nz/2.]
        ngridfile = np.zeros(ncpu + nboundary, nlevelmax)
        ngridlevel = np.zeros(ncpu + nlevelmax)

        if (nboundary > 0):
            ngridbound = np.zeros(nboundary, nlevelmax)
        # info file is needed. (and is ready)
        # cpu_list & hilbert_key : sim.info.

        try:
            lmax
        except:
            lmax = nlevelmax
        print(' >>> working resolution (lmax) =', lmax)

        xxmin = ranges[0][0]
        xxmax = ranges[0][1]
        yymin = ranges[1][0]
        yymax = ranges[1][1]
        zzmin = ranges[2][0]
        zzmax = ranges[2][1]

        dmax = max([xxmax-xxmin, yymax-yymin, zzmax-zzmin])
        for ilevel in range(1, lmax):
            dx = 0.5**ilevel
            if (dx <= dmax):
                lmin = ilevel
                break

        bit_length = lmin - 1
        maxdom = 2**bit_length
        imin = 0
        imax = 0
        jmin = 0
        jmax = 0
        kmin = 0
        kmax = 0

        if (bit_length >= 0):
            imin = int(xxmin * maxdom)
            imax = imin + 1
            jmin = int(yymin * maxdom)
            jmax = jmin + 1
            kmin = int(zzmin * maxdom)
            kmax = kmin + 1

        dkey = (2**(nlevelmax+1)/maxdom)**(ndim)
        ndom = 1
        if (bit_length > 0):
            ndom = 8
        idom = np.zeros(9)
        jdom = np.zeros(9)
        kdom = np.zeros(9)
        idom[1] = imin; idom[2] = imax; idom[3] = imin; idom[4] = imax
        idom[5] = imin; idom[6] = imax; idom[7] = imin; idom[8] = imax
        jdom[1] = jmin; jdom[2] = jmin; jdom[3] = jmax; jdom[4] = jmax
        jdom[5] = jmin; jdom[6] = jmin; jdom[7] = jmax; jdom[8] = jmax
        kdom[1] = kmin; kdom[2] = kmin; kdom[3] = kmin; kdom[4] = kmin
        kdom[5] = kmax; kdom[6] = kmax; kdom[7] = kmax; kdom[8] = kmax

        bounding_min = np.zeros(9)
        bounding_max = np.zeros(9)

        for i in range(1, ndom + 1):
            if bit_length > 0:
                order_min = self._hilbert3d([idom[i]], [jdom[i]], [kdom[i]],
                                            bit_length, 1)
                # order_min, array or single variable??
                # Will it be vectorized??
            else:
                order_min = 0

            order_min = np.asarray(order_min)
            bounding_min[i] = [order_min][0] * dkey
            bounding_max[i] = ([order_min][0] + 1) * dkey
            # [(x + 1) * dkey for x in order_min]

        cpu_min = np.zeros(9, dtype=np.int)
        cpu_max = np.zeros(9, dtype=np.int)
        bound_key = info.hilbertkey[0]

        cpu_list = np.zeros(ncpu+1, dtype=np.int)

        for impi in range(1, ncpu + 1):
            for i in range(1, ndom + 1):
                if ((bound_key[impi - 1] <= bounding_min[i]) and
                   (bound_key[impi] > bounding_min[i])):
                    cpu_min[i] = impi
                if ((bound_key[impi - 1] < bounding_max[i]) and
                   (bound_key[impi] >= bounding_max[i])):
                    cpu_max[i] = impi

        ncpu_read = 0
        cpu_read = np.zeros(ncpu + 1, dtype=np.int)
        for i in range(1, ndom + 1):
            for j in range(cpu_min[i], cpu_max[i] + 1):
                # np.arange(10,10) = [], np.arange(10,11) = [10]
                if cpu_read[j] == 0:
                    ncpu_read += 1
                    cpu_list[ncpu_read] = j
                    cpu_read[j] = 1

        for ilevel in range(1, lmax + 1):
            nx_full = 2**ilevel
            ny_full = 2**ilevel
            nz_full = 2**ilevel
            imin = xxmin * (nx_full) + 1
            imax = xxmax * (nx_full) + 1
            jmin = yymin * (ny_full) + 1
            jmax = yymax * (ny_full) + 1
            kmin = zzmin * (nz_full) + 1
            kmax = zzmax * (nz_full) + 1

# Sort cpu_list in descending npart order for memory efficiency
        cpu_list = cpu_list[cpu_list > 0]  # crop empty part

        return np.sort(cpu_list)

    def _hilbert3d(self, x, y, z, bit_length, npoint):
        '''
        Calculate hilbert doamin decomposition
        '''

        state_diagram = np.zeros((8, 2, 12), dtype=np.int)
        state_diagram[:, 0, 0] = [1, 2, 3, 2, 4, 5, 3, 5]
        state_diagram[:, 1, 0] = [0, 1, 3, 2, 7, 6, 4, 5]
        state_diagram[:, 0, 1] = [2, 6, 0, 7, 8, 8, 0, 7]
        state_diagram[:, 1, 1] = [0, 7, 1, 6, 3, 4, 2, 5]
        state_diagram[:, 0, 2] = [0, 9, 10, 9, 1, 1, 11, 11]
        state_diagram[:, 1, 2] = [0, 3, 7, 4, 1, 2, 6, 5]
        state_diagram[:, 0, 3] = [6, 0, 6, 11, 9, 0, 9, 8]
        state_diagram[:, 1, 3] = [2, 3, 1, 0, 5, 4, 6, 7]
        state_diagram[:, 0, 4] = [11, 11, 0, 7, 5, 9, 0, 7]
        state_diagram[:, 1, 4] = [4, 3, 5, 2, 7, 0, 6, 1]
        state_diagram[:, 0, 5] = [4, 4, 8, 8, 0, 6, 10, 6]
        state_diagram[:, 1, 5] = [6, 5, 1, 2, 7, 4, 0, 3]
        state_diagram[:, 0, 6] = [5, 7, 5, 3, 1, 1, 11, 11]
        state_diagram[:, 1, 6] = [4, 7, 3, 0, 5, 6, 2, 1]
        state_diagram[:, 0, 7] = [6, 1, 6, 10, 9, 4, 9, 10]
        state_diagram[:, 1, 7] = [6, 7, 5, 4, 1, 0, 2, 3]
        state_diagram[:, 0, 8] = [10, 3, 1, 1, 10, 3, 5, 9]
        state_diagram[:, 1, 8] = [2, 5, 3, 4, 1, 6, 0, 7]
        state_diagram[:, 0, 9] = [4, 4, 8, 8, 2, 7, 2, 3]
        state_diagram[:, 1, 9] = [2, 1, 5, 6, 3, 0, 4, 7]
        state_diagram[:, 0, 10] = [7, 2, 11, 2, 7, 5, 8, 5]
        state_diagram[:, 1, 10] = [4, 5, 7, 6, 3, 2, 0, 1]
        state_diagram[:, 0, 11] = [10, 3, 2, 6, 10, 3, 4, 4]
        state_diagram[:, 1, 11] = [6, 1, 7, 0, 5, 2, 4, 3]

        i_bit_mask = np.zeros(3 * bit_length)
        ind = np.arange(bit_length)
#       order      = dblarr(npoint)

        for ip in range(npoint):  # check if loop range is valid.
            # convert to binary
            x_bit_mask = self._btest(x[ip], bit_length - 1, True)
            y_bit_mask = self._btest(y[ip], bit_length - 1, True)
            z_bit_mask = self._btest(z[ip], bit_length - 1, True)

            # interleave bits
            i_bit_mask[3 * ind + 2] = x_bit_mask
            i_bit_mask[3 * ind + 1] = y_bit_mask
            i_bit_mask[3 * ind + 0] = z_bit_mask

            # build Hilbert ordering using state diagram
            cstate = 0
            # from bit_length -1 to 0 in descending order.
            for i in range(bit_length - 1, -1, -1):
                b2 = 0
                if (i_bit_mask[3 * i + 2]):
                    b2 = 1
                b1 = 0
                if (i_bit_mask[3 * i + 1]):
                    b1 = 1
                b0 = 0
                if (i_bit_mask[3 * i + 0]):
                    b0 = 1
                sdigit = b2 * 4 + b1 * 2 + b0

                nstate = state_diagram[sdigit, 0, cstate]
                hdigit = state_diagram[sdigit, 1, cstate]
                i_bit_mask[3 * i + 2] = self._btest(hdigit, 2, p_all=False)
                i_bit_mask[3 * i + 1] = self._btest(hdigit, 1, p_all=False)
                i_bit_mask[3 * i + 0] = self._btest(hdigit, 0, p_all=False)
                cstate = nstate

            # save Hilbert key as double precision real
            order = [0.0] * npoint
            for i in range(3 * bit_length):
                b0 = 0
                if (i_bit_mask[i]):
                    b0 = 1
                order[ip] += b0 * (2**i)
        return order

    def _btest(self, tmp, bit, p_all=False):
        nbit = bit
        if (not p_all and tmp != 0):
            tmp2 = int(np.log2(tmp))+1
            if (tmp2 > nbit):
                nbit = tmp2

        res = [0]*(nbit+1)

        for j in np.arange(nbit + 1, 0, -1) - 1:
            res[j] = np.int(tmp / 2**j)
            tmp -= res[j] * 2**j

        if (p_all):
            return(res)
        else:
            return(res[bit])

#%%

# base = '/home/hoseung/Work/data/AGN2/'
#base = '/home/hoseung/Work/data/Lmax_136/'
# add a method to automatically set ranges
#xrs = [[0.2, 0.23], [0.1, 0.2], [0.15, 0.23]]

#s = Sim(136, base)  #, ranges=xrs)
# s = sim.load.Sim(132, base + 'snapshots/')
# s.set_ranges(xrs)

#s.add_hydro()
#s.hydro.amr2cell(verbose=False)