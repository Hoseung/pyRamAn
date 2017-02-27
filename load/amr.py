# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:44:27 2015

@author: hoseung
"""
import numpy as np
from load.utils import read_header, read_fortran, skip_fortran

class AmrHeader():
    """
    Find the coefficients of a polynomial with the given sequence of roots.

    Returns the coefficients of the polynomial whose leading coefficient
    is one for the given sequence of zeros (multiple roots must be included
    in the sequence as many times as their multiplicity; see Examples).
    A square matrix (or array, which will be treated as a matrix) can also
    be given, in which case the coefficients of the characteristic polynomial
    of the matrix are returned.

    Parameters
    ----------
    seq_of_zeros : array_like, shape (N,) or (N, N)
        A sequence of polynomial roots, or a square array or matrix object.

    Returns
    -------
    c : ndarray
        1D array of polynomial coefficients from highest to lowest degree:

        ``c[0] * x**(N) + c[1] * x**(N-1) + ... + c[N-1] * x + c[N]``
        where c[0] always equals 1.

    Raises
    ------
    ValueError
        If input is the wrong shape (the input must be a 1-D or square
        2-D array).

    See Also
	"""
    def __init__(self):
        pass

    def _read_amr_header(self, f, skip_header = False):
        """
        Make this visible from outside, and more general
        this can be used everytime you need to skip header
        where value assigning is unnecessary
        recieve file object f rather than opening one internally.
        Are all header entries global quantaties?
        or varies with cpus?
        Local values must be separated from global values.

        parameters
        ----------
        skip_header :
            makes the core more tolerant.
            AMR header structure may change depending on the
            type of the simulation.
        """
        h1 = read_header(f, np.dtype(
                             [('ncpu', 'i4'), ('ndim', 'i4'),
                              ('ng', 'i4', (3,)), ('nlevelmax', 'i4'),
                              ('ngridmax', 'i4'), ('nboundary', 'i4'),
                              ('ngrid', 'i4'), ('boxlen', 'f8'),
                              ('outputs', 'i4', (3,))]))  # Global

        self.ncpu = h1['ncpu']
        self.ndim = h1['ndim']
        self.ng = h1['ng']
        self.nlevelmax = h1['nlevelmax']
        self.ngridmax = h1['ngridmax']
        self.ngridtot = 0
        self.nboundary = h1['nboundary']
        self.ngrid = h1['ngrid']
        self.boxlen = h1['boxlen']
        self.nnouts = h1['outputs'][0]
        self.iout = h1['outputs'][1]
        self.ifout = h1['outputs'][2]
        # Basic information that are required to read the header further.
        ncpu = h1['ncpu']
        nnouts = h1['outputs'][0]
        nlevelmax = h1['nlevelmax']
        ncoarse = np.product(h1['ng'][:])  # multiply all elements in an array

        # Global
        dtype_h2 =  np.dtype(
                     [('tout', 'f8', (nnouts,)),
                      ('aout', 'f8', (nnouts,)),
                      ('t', 'f8'),
                      ('dtold', 'f8', (nlevelmax,)),
                      ('dtnew', 'f8', (nlevelmax,)),
                      ('nsteps', 'i4', (2,))])

        dtype_h20 = np.dtype([('cmr', 'f8', (3,)),
                              ('omlkbhal', 'f8', (7,)),
                              ('expepot', 'f8', (5,)),
                              ('mass_sph', 'f8'),
                              ('headl', 'i4', (nlevelmax, ncpu,)),
                              ('taill', 'i4', (nlevelmax, ncpu,)),
                              ('numbl', 'i4', (nlevelmax, ncpu,))])

        if skip_header:
            return
            for i in range(len(dtype_h2)):
               skip_fortran(f)
            h2 = np.empty(1,dtype_h2)

            for i in range(len(dtype_h20)):
               skip_fortran(f)
            h20 = np.empty(1,dtype_h20)
            skip_fortran(f)
        else:
            h2 = read_header(f, dtype_h2)
            h20= read_header(f, dtype_h20)
            self.numbl = h20["numbl"]
            skip_fortran(f)

        if (h1['nboundary'] > 0):
            h3 = read_header(f, np.dtype([('headb', 'i4', (nlevelmax, ncpu,)),
                                          ('tailb', 'i4', (nlevelmax, ncpu,)),
                                          ('numbb', np.int16, (nlevelmax, ncpu,))]))
            self.headb = h3['headb']
            self.tailb = h3['tailb']
            self.numbb = h3['numbb']

        h4 = read_header(f, np.dtype([('htnm1m2', 'i4', (5,)),
                                      ('ordering', 'a128',(1,))]),check=True)
        self.headf = h4['htnm1m2'][0]
        self.tailf = h4['htnm1m2'][1]
        self.numbf = h4['htnm1m2'][2]
        self.used_mem = h4['htnm1m2'][3]
        self.used_mem_tot = h4['htnm1m2'][4]
        self.ordering = h4['ordering']


    # When reading 2D array, beware that fortran file is written in
    # ???-major order but python will save it in ???-major order

    # For example,
    # numbl is (ncpu x nlevelmax) array in fortran
    # and is accessed by numbl[icpu,ilevel]

    # where is the 170 coming from?


        #h4 = read_header(f, np.dtype([('bound_key', 'f8', (ncpu+1,))]),check=False)
        # Get the data type by calculating precision from the fortran block header
        alen = np.fromfile(f, np.dtype('i4'), 1)
        if alen/(ncpu + 1) == 8:
            dtype = np.float64
        elif alen/(ncpu + 1) == 16:
            dtype = np.float128
        else:
            raise Exception('Failed to detect bound_key precision.')
        self.bound_key = np.fromfile(f, dtype, ncpu + 1)
        np.fromfile(f, np.dtype('i4'), 1) # skip tail.


        h4 = read_header(f, np.dtype([('son', 'i4', (ncoarse,)),
                                      ('flag1', 'i4', (ncoarse,)),
                                      ('cpu_map', 'i4', (ncoarse,))]),check=True)

# Aquarius data has 16Byte "bound_key".
# Because of QUADHILBERT??

# check=False => Even if user gives a wrong size,
# it still reads based on what fortran binary says.

        # if assign

        self.tout = h2['tout']
        self.aout = h2['aout']
        self.t = h2['t']
        self.dtold = h2['dtold']
        self.dtnew = h2['dtnew']
        self.nstep = h2['nsteps'][0]
        self.nstep_coarse = h2['nsteps'][1]

        self.const = h20['cmr'][0]
        self.mass_tot0 = h20['cmr'][1]
        self.rho_tot = h20['cmr'][2]
        self.Om = h20['omlkbhal'][0]
        self.Ol = h20['omlkbhal'][1]
        self.Ok = h20['omlkbhal'][2]
        self.Ob = h20['omlkbhal'][3]
        self.h0 = h20['omlkbhal'][4]
        self.aexp_ini = h20['omlkbhal'][5]
        self.boxlen = h20['omlkbhal'][6]

        self.aexp = h20['expepot'][0]
        self.hexp = h20['expepot'][1]
        self.aexp_old = h20['expepot'][2]
        self.epot_tot_ini = h20['expepot'][3]
        self.epot_tot_old = h20['expepot'][4]

        self.mass_sph = h20['mass_sph']

        self.headl = h20['headl']
        self.taill = h20['taill']
        self.numbl = h20['numbl']
#        self.numbot = h2['numbot'] # This value has been skipped

        self.son = h4['son']
        self.flag1 = h4['flag1']
        self.cpu_map = h4['cpu_map']

class Grid():
    def __init__(self):
        self.ncpu = 0
        self.ndim = 0
        self.time = 0.
        self.aexp = 0.
        self.nlevelmax=0
        self.boxlen = 0.
        self.ngridtot=0
        self.ngridarr=None
        self.levellist=None

    def set_data(self, ncpu=0, ndim=0, time=0., aexp=0.,
                nlevelmax=0, boxlen=0., ngridtot=0,
                ngridarr=None, levellist=None):
        self.ncpu = ncpu
        self.ndim = ndim
        self.time = time
        self.aexp = aexp
        self.nlevelmax=nlevelmax
        self.boxlen = boxlen
        self.ngridtot=ngridtot
        self.ngridarr=ngridarr
        self.levellist=levellist


class Amr():
    """
    AMR class, which is required by Hydro class.

    """

    def __init__(self, info, cpus=None, load=True):
        import os
        """
        Parameters
        ----------
        info : load.info.Info class

        """
        snout = str(info.nout).zfill(5)
        self.info = info
        self.cpus = cpus

        self._fnbase = os.path.join(info.base, info.data_dir) + 'output_' + snout + '/amr_' + snout + '.out'
        try:
            f = open(self._fnbase + '00001', "rb")
        except:
            import glob
            amrs = glob.glob(self._fnbase + "*")
            f = open(amrs[0], "rb")

        self.header = AmrHeader()
        self.header._read_amr_header(f)
        if load:
            self.load()
        f.close()

    def _load_mesh(f, ndim=3):
        for i in np.arange(ndim):
            read_fortran(f, np.dtype('f8'))

    def get_zoomin(self, scale=1.0):
        """
        Returns a spherical region encompassing maximally refined cells.

        Parameters
        ----------
        scale : float
            The radius of the returned sphere is scaled by 'scale'.

        """
        from utils import sampling
        imin = np.where(self.dm['m'] == self.dm['m'].min())
        xr = [self.dm['px'][imin].min(), self.dm['px'][imin].max()]
        yr = [self.dm['py'][imin].min(), self.dm['py'][imin].max()]
        zr = [self.dm['pz'][imin].min(), self.dm['pz'][imin].max()]

        xc = 0.5 * sum(xr)
        yc = 0.5 * sum(yr)
        zc = 0.5 * sum(zr)

        radius = 0.5 * max([xr[1]-xr[0], yr[1]-yr[0], zr[1]-zr[0]]) * scale
        #print(radius)
        return sampling.set_region(centers=[xc, yc, zc], radius=radius)


    def load(self, verbose=False):
        """
        Load AMR data.

        Need info file.

        Notes
        -----
        The building block of FTT AMR is an oct, which is a group of 8 cells and data
        either associate with cells or the oct(called 'mesh' in IDL analysis routine).
        An Oct consists of (level, xyz coordinates, pointer to the parent,
        pointer to 6 neighbouring parents, pointer to the child oct)

        cpu map and refinement map is additionaly needed to restart a simulation.
        octs of the same level are written at once. (So you need to loop over ilevel)

        <Todo>
        Think about when you will need to load the amr file.
        It's rather unclear.

        -> When I want to draw cell map!
        """

        # global header variables are available.
        icpu = 0  # icpu
        cpus = self.cpus
        ndim = self.header.ndim

        ncpu = self.header.ncpu
        nboundary = self.header.nboundary
        nlevelmax = self.header.nlevelmax

        ncell = 2**ndim
        xbound=np.zeros(3)
        # Size of arrays

        if icpu == 0:
            listmax = ncpu
        elif icpu > 0:
            listmax = ncpu + nboundary

        ngridarr = np.zeros((nlevelmax, listmax), dtype=np.int32)
        #ngridarr = np.zeros((nlevelmax, listmax))
        levellist = [[0] * listmax for i in range(nlevelmax)] # nlevelmax by listmax list.
        #levellist = [[0] * nlevelmax for i in range(listmax)] # nlevelmax by listmax list.
        #np.zeros((nlevelmax, listmax), dtype=np.int32)
        llist = 0
        self.header.ngridtot = 0
        for jcpu in cpus:
            if(verbose):
                self.print_cpu(self, icpu)

            #print(self._fnbase + str(jcpu).zfill(5))
            f = open(self._fnbase + str(jcpu).zfill(5), "rb")  # +1

            # read header
            header_icpu = AmrHeader()
            header_icpu._read_amr_header(f)
            #print(header_icpu.ngridtot)
            #print(header_icpu.ngrid)
            self.header.ngridtot += header_icpu.ngrid

            numbl = header_icpu.numbl
            if nboundary > 0:
                numbb = header_icpu.numbb  # None if nboundary = 0
                xbound[:] = nx/2., ny/2., nz/2.
            else:
                numbb = np.zeros(np.shape(numbl)) # need an empty array.

            ngridtot = 0

            kcpumin = 1
            kcpumax = nboundary + ncpu
            nlevel = 0

            for ilevel in np.arange(nlevelmax):
                for kcpu in np.arange(kcpumin, kcpumax + 1):

                    if (kcpu <= ncpu):
                        ng = numbl[ilevel][kcpu-1]
                    else:
                        ng = numbb[ilevel][kcpu - ncpu - 1]
                    if icpu == 0:
                        if kcpu == jcpu:
                            ngridarr[ilevel][kcpu-1] = ng
                    else:
                        ngridarr[ilevel][kcpu-1] = ng

                    if (ng > 0):
                        ngridtot = ngridtot + ng
                        if (verbose):
                            print("Level %2d has %6d grids in proc %4d"
                                  % (ilevel + 1, ng, kcpu))
                        #print(f.tell())
                        nlevel = nlevel + 1  # number of valid (ng > 0) levels
                        mesh = {"ilevel":ilevel,
                                "nc":ng,
                                "xg":np.zeros((ndim, ng), dtype=np.float64),
                                "son":np.zeros((ng, ncell), dtype=np.int32)}

                        # Read actual data
                        i_skip = 0
                        ind_current = read_fortran(f, np.dtype('i4'), ng)
                        # row-major
                        ind_next = read_fortran(f, np.dtype('i4'), ng)
                        ind_prev = read_fortran(f, np.dtype('i4'), ng)
                        for idim in range(ndim): # gird center
                            xx = read_fortran(f, np.dtype('f8'), ng)
                            mesh["xg"][idim] = xx - xbound[idim]
                        # father index
                        read_fortran(f, np.dtype('i4'), ng)

                        for idim in range(2*ndim):  # neighbour index
                            read_fortran(f, np.dtype('i4'), ng)

                        for idim in range(2**ndim):  # son index
                            mesh["son"] = read_fortran(f, np.dtype('i4'), ng)

                        for idim in range(2**ndim):  # cpu map
                            read_fortran(f, np.dtype('i4'), ng)

                        for idim in range(2**ndim):  # ref map
                            read_fortran(f, np.dtype('i4'), ng)

                        if icpu == 0:
                            if (kcpu == jcpu):
                                levellist[ilevel][jcpu-1] = mesh
                        else:
                            levellist[ilevel][kcpu-1] = mesh

                        llist += 1
            f.close

#        self.grid = Grid()
        #self.grid.set_data(ncpu=listmax, ndim=ndim, time=t, aexp=aexp,
        self.ngridarr = ngridarr
        self.levellist = levellist
        #grid.set_data(ngrid=ngridarr, levellist=levellist)

        return
