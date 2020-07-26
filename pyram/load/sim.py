# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:44:27 2015

@author: hoseung
"""
import numpy as np
from ..utils.io_module import read_header, read_fortran, skip_fortran
from ..utils.cosmology import Timeconvert

class Simbase():
    """
    Common class to [amr, hydro, part] classes.

    """
    def __init__(self,
                 cosmo=True,
                 verbose=False):
        self.cosmo=cosmo
        self.ranges=None
        self.region=None
        self.nout=None
        self.cpu_fixed=False
        self.cpus=None
        self.verbose=verbose

    def set_ranges(self, ranges=[[0, 1], [0, 1], [0, 1]]):
        if ranges is not None:
            if not hasattr(self, "amr"):
                self.add_amr(load=False)

            nr = np.asarray(ranges) # Now it is a class
            if not(nr.shape[0] == 3 and nr.shape[1] == 2):
                # Because actual operation on the given input(ranges)
                # does not take place soon, it's not a good place to use
                # try & except clause. There is nothing to try yet.
                print(' Error!')
                print('Shape of ranges is wrong:', nr.shape)
                print('example : [[0.1,0.3],[0.2,0.4],[0.6,0.8]] \n')
            else:
                self.ranges = ranges
                self.set_cpus()
        else:
            self.ranges = None

    def add_amr(self, load=False):
        self.amr = Amr(info=self.info, cpus=self.cpus, load=load)
        if self.verbose: 
            print("An AMR instance is created\n", self.__class__.__name__)

    def set_cpus(self, cpus=None, lock=False):
        """
        Determine the list of CPUs for the range.
        Requires self.range.
        """
        if cpus is None:
            cpus = self._hilbert_cpulist(self.info, self.ranges,
                                         amrheader=self.amr.header)
        if not self.cpu_fixed:
            self.cpus = np.array(cpus)
        # Lock, but only if there is a valid list of cpus.
        if lock and self.cpus is not None:
            self.cpu_fixed = True

    def unlock_cpus(self):
        self.cpu_fixed=False

    def show_cpus(self):
        print(" ncpus : %s \n" % self.cpus)

    def _hilbert_cpulist(self, info, ranges, amrheader=None):
        """
        After determining cpus, read the cpu files and cut off data
        that are outside ranges.
        -> cpu files contain data points within a given ranges
        BUT! they also contain other data points outside the ranges.
        """
        from .a2c import a2c

        if amrheader is None:
            if not(hasattr(self, 'amr')):
                print("[sim._hilbert_cpulist] No AMR instance,")
                print("[sim._hilbert_cpulist] Loading one...")
                self.add_amr(load=False)# load only meta data
            amrheader = self.amr.header

        nlevelmax = amrheader.nlevelmax
        #nboundary = amrheader.nboundary
        ndim = self.info.ndim
        ncpu = self.info.ncpu_tot  #ncpu_tot
        # header of part output is needed (and is prepared.)

        # header of amr output is needed.
        #two2ndim = 2**ndim
        #nx = amrheader.ng[0]
        #ny = amrheader.ng[1]
        #nz = amrheader.ng[2]

        #xbound = [nx/2., ny/2., nz/2.]
        #ngridfile = np.zeros(ncpu + nboundary, nlevelmax)
        #ngridlevel = np.zeros(ncpu + nlevelmax)

        #if (nboundary > 0):
        #    ngridbound = np.zeros(nboundary, nlevelmax)

        try:
            lmax
        except:
            lmax = nlevelmax

        xxmin = ranges[0][0]
        xxmax = ranges[0][1]
        yymin = ranges[1][0]
        yymax = ranges[1][1]
        zzmin = ranges[2][0]
        zzmax = ranges[2][1]

        dmax = max([xxmax-xxmin, yymax-yymin, zzmax-zzmin])
        # lmin = lmax # sometimes even smallest dx is larger than the region size.
        lmin = info.lmin # default value. But, which value should I pick?
        for ilevel in range(1, lmax):
            dx = 0.5**ilevel
            if (dx < dmax):
                lmin = ilevel
                break

        bit_length = lmin - 1
        maxdom = 2**bit_length
        imin, imax, jmin, jmax, kmin, kmax = 0,0,0,0,0,0

        if (bit_length > 0):
            imin = int(xxmin * maxdom)
            imax = imin + 1
            jmin = int(yymin * maxdom)
            jmax = jmin + 1
            kmin = int(zzmin * maxdom)
            kmax = kmin + 1

        dkey = (2**(nlevelmax+1)/maxdom)**(ndim)
        ndom = 1
        if (bit_length > 0): ndom = 8
        idom = np.zeros(9)
        jdom = np.zeros(9)
        kdom = np.zeros(9)
        idom[1] = imin; idom[2] = imax
        idom[3] = imin; idom[4] = imax
        idom[5] = imin; idom[6] = imax
        idom[7] = imin; idom[8] = imax
        jdom[1] = jmin; jdom[2] = jmin
        jdom[3] = jmax; jdom[4] = jmax
        jdom[5] = jmin; jdom[6] = jmin
        jdom[7] = jmax; jdom[8] = jmax
        kdom[1] = kmin; kdom[2] = kmin
        kdom[3] = kmin; kdom[4] = kmin
        kdom[5] = kmax; kdom[6] = kmax
        kdom[7] = kmax; kdom[8] = kmax

        bounding_min = np.zeros(9)
        bounding_max = np.zeros(9)

        for i in range(1, ndom + 1):
            if bit_length > 0:
                order_min = a2c.hilbert3d([idom[i]], [jdom[i]], [kdom[i]],
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
                    (bound_key[impi    ] >  bounding_min[i])):
                    cpu_min[i] = impi
                if ((bound_key[impi - 1] <  bounding_max[i]) and
                    (bound_key[impi    ] >= bounding_max[i])):
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

        # Sort cpu_list in descending npart order for memory efficiency
        cpu_list = cpu_list[cpu_list > 0]  # crop empty part

        return np.sort(cpu_list)

def empty_config():
    conf = {'cosmo':None,
            'sim_type':None,
            'dmo':None,
            'longint':None,
            'sink_path':'./',
            }
    return conf
"""
Config templates are in config.py. But what is actually used need to belong to a sim instance.
"""


class Sim(Simbase):
    """
    Defines the 'host class' of part, amr, hydro, info.


    Methods
    -------
    setup(self, nout=None, base='./', data_dir='snapshots/', ranges=[[0.0,1.0],[0.0,1.0],[0.0,1.0]], dmo=False)
        Setup basic parameters need for an AMR instance.

    add_info(self, load=False)

    add_hydro(self, load=False, lmax=None)

    add_part(self, ptypes=[], load=False, fortran=True, dmo=False, **kwargs)

    Examples
    --------
    >>> import load
    >>> s = load.sim.Sim(187, ranges=[[0.3,0.4], [0.35,0.45], [0.1,0.2]])

    Notes
    -----
    Global information is stored in this class:
        ndim, ncpu, base, type of simulation (DMO, zoom, and so on)
    Later it will also include .nml information.
    (Romain's git version generates such output in text files)

    Currently this class deals with single snapshot.
    But I hope to expand it for multiple snapshots.
    """
    def __init__(self, nout, base='./', data_dir='snapshots/',
                 ranges=[[0.0,1.0],[0.0,1.0],[0.0,1.0]], setup=True, 
                 dmo=False, region=None, cosmo=True, sim_type='none'):
        """
            Parameters
            ----------
            nout: int
            base: str, optional
            data_dir: str, optional
            ranges: list, optional
                [[xmin,xmax],[ymiin,ymax],[zmin,zmax]]
            dmo: bool, optional
            setup: bool, optional
            region: (region) class, optional

        """
        super(Sim,self).__init__()
        # should call parent class' init.

        self.nout = nout
        self.base = base
        # info appreciates nout and base (not mandatary, though)
        self.config = empty_config()
        self.config['dmo'] = dmo
        self.config['cosmo'] = cosmo
        self.config['sim_type'] = sim_type.lower()
        self.data_dir = data_dir
        self.add_info()
        # set_data_dir and set_range needs info instance be exist.
        self.set_ranges(ranges)
        if region is not None:
            print("setting ranges")
            ranges = [[region.xc - region.radius, region.xc + region.radius],
                      [region.yc - region.radius, region.yc + region.radius],
                      [region.zc - region.radius, region.zc + region.radius]]
            self.region = region
            print(self.region)
        else:
            self.region=None

        if setup:
            self.setup(ranges, dmo)

    def _all_set(self):
        return (self.nout is not None) & (self.base is not None)

    def setup(self, ranges=[[0.0,1.0],[0.0,1.0],[0.0,1.0]], dmo=False):
        if self.nout is None:
            raise ValueError("Note that 'nout' is not set. \n use sim.Sim.set_nout(nout)")
        self.add_info()
        # set_data_dir and set_range needs an info instance.
        self.set_ranges(ranges)

        if self._all_set():
            #self.add_amr()
            if self.ranges is not None:
                self.set_cpus(self._hilbert_cpulist(self.info, self.ranges))
        print('Simulation set up.')

    """
    Programming note.
    No need to use .gettter and .setter always.
    Just make values public if they need to be accessed.
    If an additional operation is needed on setting a variable,
    use @property macro like this.

    Here, _base is encapsulated.
    """
    @property
    def base(self):
        return self._base
    @base.setter
    def base(self, base):
        from os import path
        self._base = path.abspath(base)

    @property
    def nout(self):
        return self._nout
    @nout.setter
    def nout(self, nout):
        """
            Sets output number.
            a list of nouts will be supported in the future, soon, I hope.
        """
        self._nout = nout

    @property
    def data_dir(self):
        return self._data_dir
    @data_dir.setter
    def data_dir(self, data_dir):
        """
        By default, simulation outputs are in simulation_base/snapshots/
        """
        from os import path
        #self.data_dir =  path.join(self.base, '', data_dir, '')
        self._data_dir = path.join('', data_dir, '')

    def show_base(self):
        print("setting the base(working) directory to :", self.base)

    def add_hydro(self, load=True, lmax=None, region=None, ranges=None,
                  cpu=False, **kwargs):
        """
        Add a hydro instance to the simulation instance.

        parameters
        ----------
        lmax : int
            maximum AMR level of hydro variable retrieved.
        load : bool
            If false, an hydro instance is added without cell data.
        """
        from . import hydro
        if region is None:
            region = self.region
        if ranges is None:
            ranges = self.ranges

        self.hydro = hydro.Hydro(info=self.info,
                                 cpus=self.cpus,
                                 cpu_fixed=self.cpu_fixed,
                                 region=region,
                                 ranges=ranges)

        if load :
            if lmax is None:
                lmax = self.info.lmax
            self.hydro.load(lmax=lmax, cpu=cpu, **kwargs)
        else:
            print("Use hydro.amr2cell() to load hydro variables")

    def add_info(self, load=False):
        from . import info
        self.info = info.Info(self.nout,
                              self.base,
                              load=load,
                              data_dir = self.data_dir,
                              cosmo = self.cosmo)
        if self.cosmo:
            self.tc = Timeconvert(self.info)
            self.info.tGyr = self.tc.time2gyr(self.info.time, zred_now = self.info.zred)

    def add_part(self, ptypes=[], load=True, fortran=True, dmo=False, **kwargs):
        """
        Add a particle instance to the simulation instance.
        Requires types of particles and particle data.
        load = True  to load actual data on creating the instance

        parameters
        ----------
        ptypes : list of particle type and information.
                ["dm id pos"] or ["dm id pos", "star mass vel"]

        """
        if dmo:
            self.dmo = True
        from . import part
        print("Types of particles you want to load are: ", ptypes)

        # To do. instead of specifying every variables,
        # make the Part object inherit common variables from the father class instance, sim.
        # use inspect module??
        self.part = part.Part(info=self.info,
                              ptypes=ptypes,
                              data_dir=self.data_dir,
                              config=self.config,
                              base=self.base,
                              cpus=self.cpus,
                              cpu_fixed=self.cpu_fixed,
                              region=self.region,
                              ranges=self.ranges,
                              **kwargs)
        print("A particle instance is created\n")

        if load:
            read_metal=False
            if "metal" in ptypes:
                read_metal=True
            self.part.load(fortran=fortran, read_metal=read_metal)
        else:
            print("Use part.load() to load particle")

class AmrHeader():
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

        # numbl is (ncpu x nlevelmax) array in fortran
        # and is accessed by numbl[icpu,ilevel]

        #h4 = read_header(f, np.dtype([('bound_key', 'f8', (ncpu+1,))]),check=False)
        # Get the data type by calculating precision from the fortran block header
        self.key_size = int(np.fromfile(f, np.dtype('i4'), 1))
        if self.key_size/(ncpu + 1) == 8:
            dtype = np.float64
        elif self.key_size/(ncpu + 1) == 16:
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
        # self.numbot = h2['numbot'] # This value has been skipped

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

class Amr(Simbase):
    """
    AMR class, which is required by Hydro class.
    """

    def __init__(self, nout=None, info=None, cpus=None, load=True):
        import os
        """
        Parameters
        ----------
        info : load.info.Info class

        """
        super(Amr, self).__init__()
        if info is not None:
            self.info = info
        else:
            assert nout is not None, "either info or nout is required"
            from .info import Info
            print("[Amr.__init__] Loading info")
            self.info = Info(nout=nout)

        snout = str(self.info.nout).zfill(5)

        # Update attributes
        if cpus is not None:
            self.cpus = cpus

        self._fnbase = os.path.join(self.info.base, self.info.data_dir) +\
                    'output_' + snout + '/amr_' + snout + '.out'
        try:
            f = open(self._fnbase + '00001', "rb")
        except:
            import glob
            amrs = glob.glob(self._fnbase + "*")
            f = open(amrs[0], "rb")

        self.header = AmrHeader()
        self.header._read_amr_header(f)
        self.get_current_lmax()

        # set_data_dir and set_range needs an info instance.
        self.set_ranges()

        #if self._all_set():
            #self.add_amr()
        #if self.ranges is not None:
        #    self.set_cpus(self._hilbert_cpulist(self.info, self.ranges,
        #                                        amrheader=self.header))

        if load:
            self.load()
        f.close()

    def set_ranges(self, ranges=[[0, 1], [0, 1], [0, 1]]):
        """
        Override set_range in the Simbase class.

        """
        if ranges is not None:

            nr = np.asarray(ranges) # Now it is a class
            if not(nr.shape[0] == 3 and nr.shape[1] == 2):
                # Because actual operation on the given input(ranges)
                # does not take place soon, it's not a good place to use
                # try & except clause. There is nothing to try yet.
                print(' Error!')
                print('Shape of ranges is wrong:', nr.shape)
                print('example : [[0.1,0.3],[0.2,0.4],[0.6,0.8]] \n')
            else:
                self.ranges = ranges
                self.set_cpus(self._hilbert_cpulist(self.info, self.ranges,
                                                    amrheader=self.header))
        else:
            self.ranges = None

    #def _load_mesh(f, ndim=3):
    #    for i in np.arange(ndim):
    #        read_fortran(f, np.dtype('f8'))
    def get_current_lmax(self):
        imax_level = np.argmax(np.sum(self.header.numbl, axis=1) == 0)
        if imax_level == 0:
            self.lmax_now = self.info.lmax
        else:
            self.lmax_now = imax_level

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
        levellist = [[0] * listmax for i in range(nlevelmax)] # nlevelmax by listmax list.
        llist = 0
        self.header.ngridtot = 0
        for jcpu in cpus:
            if(verbose):
                self.print_cpu(self, icpu)

            f = open(self._fnbase + str(jcpu).zfill(5), "rb")  # +1

            # read header
            header_icpu = AmrHeader()
            header_icpu._read_amr_header(f)
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

        self.ngridarr = ngridarr
        self.levellist = levellist

        return

