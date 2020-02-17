# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 23:27:16 2015

@author: hoseung
"""

import numpy as np
from . import Simbase


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
            region: (region) class, optional

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
            ranges = [[region.xc - region.radius, region.xc + region.radius],
                      [region.yc - region.radius, region.yc + region.radius],
                      [region.zc - region.radius, region.zc + region.radius]]
        else:
            self.region=None

        if setup:
            self.setup(nout, base, ranges, dmo)

    def _all_set(self):
        return (self.nout is not None) & (self.base is not None)

    def setup(self, nout=None, base='./',
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
        from . import info
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
        from . import part
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
