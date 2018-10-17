# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 23:27:16 2015

@author: hoseung
"""

import numpy as np

class Simbase():
    """
    base
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

    def add_amr(self, load=False):
        from load.amr import Amr
        self.amr = Amr(self.info, cpus=self.cpus, load=load)
        if self.verbose: print("An AMR instance is created\n", self.__class__.__name__)

    def unlock_cpus(self):
        self.cpu_fixed=False

    def set_cpus(self, cpus=None, lock=False):
        """
        Determine the list of CPUs for the range.
        Requires self.range.
        """
        if cpus is None:
            cpus = self._hilbert_cpulist(self.info, self.ranges)

        if not self.cpu_fixed:
            self.cpus = np.array(cpus)
        # Lock, but only there is a valid list of cpus.
        if lock and self.cpus is not None:
            self.cpu_fixed = True

    def show_cpus(self):
        print(" ncpus : %s \n" % self.cpus)

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


    def _hilbert_cpulist(self, info, ranges):
        """
        After determining cpus, read the cpu files and cut off data
        that are outside ranges.
        -> cpu files contain data points within a given ranges
        BUT! they also contain other data points outside the ranges.
        """
        from load.a2c import a2c
        if not(hasattr(self, 'amr')):
            print("[sim._hilbert_cpulist] No AMR instance,")
            print("[sim._hilbert_cpulist] Loading one...")
            self.add_amr(load=False)# load only meta data

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
#        lmin = lmax # sometimes even smallest dx is larger than the region size.
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
                 ranges=[[0.0,1.0],[0.0,1.0],[0.0,1.0]], dmo=False,
                 setup=True, region=None, cosmo=True):
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
            region: (region)dict, optional

        """
        super(Sim,self).__init__()
        # should call parent class' init.

        self.nout = nout
        self.base = base
        # info appreciates nout and base (not mandatary, though)
        self.dmo = dmo
        self.cosmo = cosmo
        self.data_dir = data_dir
        self.add_info()
        # set_data_dir and set_range needs info instance be exist.
        self.set_ranges(ranges)
        if region is not None:
            ranges = [[region['xc'] - region['radius'], region['xc'] + region['radius']],
                      [region['yc'] - region['radius'], region['yc'] + region['radius']],
                      [region['zc'] - region['radius'], region['zc'] + region['radius']]]
        else:
            self.region=None

        if setup:
            self.setup(nout, base, data_dir, ranges, dmo)

    def _all_set(self):
        return (self.nout is not None) & (self.base is not None)

    def setup(self, nout=None, base='./', data_dir='snapshots/',
                 ranges=[[0.0,1.0],[0.0,1.0],[0.0,1.0]], dmo=False):
        self.nout = nout
        self.base = base
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
        from load import hydro
        if region is None:
            region = self.region
        if ranges is None:
            ranges = self.ranges

        self.hydro = hydro.Hydro(info=self.info,
                                 cpus=self.cpus,
                                 cpu_fixed=self.cpu_fixed,
                                 region=region,
                                 ranges=ranges)
        #for kwarg in kwargs:
        #    if kwarg== "amr2cell_params":
        #        amr2cell_params = kwarg
        #print(amr2cell_params)
        if load :
            if lmax is None:
                lmax = self.info.lmax
            self.hydro.load(lmax=lmax, cpu=cpu, **kwargs)
        else:
            print("Use hydro.amr2cell() to load hydro variables")

    def add_info(self, load=False):
        from load import info
        self.info = info.Info(self.nout,
                              self.base,
                              load=load,
                              data_dir = self.data_dir,
                              cosmo = self.cosmo)
#        self.info.setup()

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
        from load import part
        print("Types of particles you want to load are: ", ptypes)

        # To do. instead of specifying every variables,
        # make the Part object inherit common variables from the father class instance, sim.
        # use inspect module??
        self.part = part.Part(info=self.info,
                              ptypes=ptypes,
                              data_dir=self.data_dir,
                              dmo=self.dmo,
                              base=self.base,
                              cpus=self.cpus,
                              cpu_fixed=self.cpu_fixed,
                              region=self.region,
                              ranges=self.ranges,
                              cosmo=self.cosmo, **kwargs)
        print("A particle instance is created\n")

        if load:
            read_metal=False
            if "metal" in ptypes:
                read_metal=True
            self.part.load(fortran=fortran, read_metal=read_metal)
        else:
            print("Use part.load() to load particle")


    def search_zoomin_region(self, *args, **kwargs):
        """
        Determine Zoomin region.

        Returns a spherical region encompassing maximally refined cells.

        Notes
        -----
        If part is not given, load particles.

        Not only part, but also hydro or amr can be used to find the zoomin region!
        But, don't want to write a new code.
        """
        if hasattr(self, 'part'):
            print("Have part")
            self.set_zregion(self.part.search_zoomin( *args, **kwargs))

        if hasattr(self, 'amr'):
            print("Have amr")
            self.set_zregion(self.amr.search_zoomin( *args, **kwargs))

    def set_zregion(self, zregion):
        """
        Set zoom-in region for Zoom-in simulations.
        """
        self.zregion = zregion
        self.info.zregion = zregion
