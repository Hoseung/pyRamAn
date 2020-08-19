# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 22:15:03 2015

@author: hoseung
"""
import numpy as np
from .sim import Simbase
from ..utils.io_module import read_header, read_fortran, skip_fortran
from .info import Info
from .part_load import part_load_module
from .readr import readr
from numpy.core.records import fromarrays as fromarrays
from ..config import part_dtype
from os.path import join, exists
from ..config import sinkprop_format, sink_prop_dtype_drag, sink_prop_dtype, sink_prop_dtype_drag_fornax


# Check if there is dm / star / sink in quantities list
# or mass, id, vel, and so on in ptype list.
# Only exact match (lower or upper case) works. No tolerence for errata.

def compute_boundary(cpumap, cpulist):
    """
    I think the name is misleading... 
    """
    bound = np.searchsorted(cpumap, cpulist)
    return np.concatenate([bound, [cpumap.size]])

class Ptype():
    def __init__(self):
        self.mass = False
        self.id = False
        self.vel = False
        self.pos = False
        self.ref = False
        self.time = False
        self.metal = False

    def _set_quantities(self, quantities):
        qlist = ["mass", "pos", "id", "vel", "ref", "time", "metal"]
        for i in qlist:
            if i in quantities:
                setattr(self, i, True)

class Ptypes():
    def __init__(self, pt, pq):
        for it, iq in zip(pt, pq):
            if it == "dm":
                self._add_dm()
                self.dm._set_quantities(iq)
                self.load_dm = True
            if it == "star":
                self._add_star()
                self.star._set_quantities(iq)
                self.load_star = True
            if it =="sink":
                self._add_sink()
                self.sink._set_quantities(iq)
                self.load_sink = True

    def _add_dm(self):
        self.dm = Ptype()

    def _add_star(self):
        self.star = Ptype()

    def _add_sink(self):
        self.sink = Ptype()

class Part(Simbase):
    """
    Particle data container supports different types of particles
    and some meta data.

    DM, star, sink
    """
    # Rule of thumb: initialize everything in __init__
    # Otherwise, a method may fail to find an object to act on.
    # This is called 'Object consistency'
    #
    # Part inherits from Simbase, but if I define __init__ once again here,
    # the __init__ methods of Simbase and Part does not merge.
    # (Think about it, that's very strange).
    # Instead, Simbase.__init__ is overridden.

    def __init__(self, config, nout=None, info=None, 
                ptypes=None, base='./',
                region=None, ranges=None, cpus=None,
                cpu_fixed=False,
                data_dir='snapshots/',
                dmref=False,
                dmvel=False,
                dmmass=True,
                load=False,
                fortran=True,
                verbose=False):
        """
        parameters
        ----------
        config : Dict with simulation type information.
            ['cosmo':[True, False], 'part':['yzics', 'nh', 'fornax', ...]
        ptypes : list of particle types and quantities.
                ["dm id pos"] or ["dm id pos", "star mass vel"]
        dmo : logical
            If True, a faster, DMO read routine invoked (NOT distingushing particle types).
        region : region dict
            only part of snapshot is loaded.
        dmref : logical
            Set True if the snapshot has DM ref information.
        dmvel : logical

        fortran : logical
            use fortran function to load data. True by Default
        verbose : logical
        Notes
        -----
        info is required for domain decomposition.

        """
        super(Part, self).__init__()
        self.config=config # {'part':'2017', 'cosmo':True}
        self.longint = False
        self.dtype = part_dtype[self.config.sim_type]
        if info == None:
            assert not nout == None, "either info or nout is required"
            info = Info(base=base,nout=nout, cosmo=self.config.cosmo, data_dir=data_dir)
        self.info = info
        self.nout = info.nout
        # To do: Need to work out cpu list mechanism
        self.cpus = cpus
        self.cpu_fixed=cpu_fixed
        self.classic_format = self.config.sim_type not in (['fornax', 'ng'])
        try:
            self.ncpu = len(self.cpus)
        except:
            self.ncpu = 0
        self.nstar = 0
        self.nsink = 0

        if verbose: print("Part", base)
        try:
            self.base = self.info.base
            if verbose: print(self.info.base)
        except:
            self.base = base
        
        #self.data_dir = data_dir
        
        self.pt = []
        self.pq = []
        self.pqset = set([])
        self.dmo = config.dmo
        self.dm_with_ref = dmref
        self.dm_with_vel = dmvel
        self.dm_with_mass = dmmass

        self.set_fnbase(self.base, self.info.data_dir, 'part')

        if region is not None:
            ranges = region.ranges
        if ranges is not None:
            self.set_ranges(ranges=ranges)
        else:
            try:
                self.set_ranges(ranges=self.info.ranges)
            except:
                pass
                # If range, reigon, info.ranges are all None,
                # then the region information is meant to be omitted.
                # probably icpu=[1,2,3] option is used.

        # header structure
        if self.config.sim_type != "fornax":
            self._ramses_particle_header = np.dtype([('ncpu', 'i4'),
                                                    ('ndim', 'i4'),
                                                    ('npart', 'i4'),
                                                    ('randseed', 'i4', (4,)),
                                                    ('nstar', 'i4'),
                                                    ('mstar', 'f8'),
                                                    ('mstar_lost', 'f8'),
                                                    ('nsink', 'i4')])
            self._get_basic_info()
            # Depending on user's choice, generate dm, star, sink classes

        if ptypes is not None: self.set_dtype(ptypes)
        if load: self.load(fortran=fortran)

    def set_dtype(self, ptypes):
        self.setwhattoread(ptypes)
        if self.classic_format:
            # check if star particle exists
            dtype_name = self.config.sim_type
            if self.nstar == 0: #NOPE, has_star_this_snapshot
                dtype_name += '_dm_only'
            self.rur_dtype = part_dtype[dtype_name]
            
        elif(self.longint):
            if(self.config.sim_type in ['iap', 'gem', 'nh']):
                self.rur_dtype = part_dtype['gem_longint']
        else:
            self.rur_dtype = part_dtype[self.config.sim_type]        

    def _read_nstar(self):
        """
        use _get_basic_info instead.
        """
        fn = self.base + f'/snapshots/output_{self.nout:05d}/part_{self.nout:05d}.out00001'
        with open(fn, 'rb') as part_file:
            part_file.skip_records(4)
            if(not self.longint):
                return part_file.read_ints()
            else:
                return part_file.read_longs()

    def mass2msun(self):
        """
        No use case?
        """
        for ptype in self.pt:
            part = getattr(self, ptype)
            if max(part["m"]) < 100:
                part["m"] *= self.info.msun

    def setwhattoread(self, ptypes):
        """
        Sorts out which types of quantities to read.
        
        Note
        ----
        Because there is no distinction b/w particles in old format, you need to read all ID
        even if you want only DM IDs.
        """

        self.family_keys = []
        for pp in ptypes:
            pp = pp.lower()
            try:
                self.family_keys.append(int(pp.split()[0]))
                pp = " ".join(pp.split()[1:])
            except:
                pass
            self.pt.append(pp.split()[0])
            self.pq.append(pp.split()[1:])
            self.pqset.update(pp.split()[1:])

        self.ptypes = Ptypes(self.pt, self.pq)
        # Check if there is dm / star / sink in quantities list
        # or mass, id, vel, and so on in ptype list.
        # Only exact match (lower or upper case) works. No tolerence for errata.
        quantities = set(["mass", "id", "vel", "ref", "time", "metal"])
        self.pqset.intersection_update(quantities)
        # any elements not listed in qunatities are dropped.

    def setDmQuantities(self, vel=False, mass=True, ref=False):
        self.dm_with_vel = vel
        self.dm_with_ref = ref

    def _get_basic_info(self):
        try:
            f = open(self._fnbase + '00001', "rb")
        except:
            from glob import glob
            parts = glob(self._fnbase + "*")
            f = open(parts[0], "rb")

        header = read_header(f, self._ramses_particle_header)
        self.ncpu = header['ncpu']
        self.nstar = header['nstar']
        self.nsink = header['nsink']

    def _get_npart_arr(self, cpus):
        npart_arr = []
        for icpu in cpus:
            with open(self._fnbase + str(icpu).zfill(5), "rb") as f:
                header = read_header(f, self._ramses_particle_header)
                npart_arr.append(header['npart'])

        return npart_arr  # ,nstar_this, my_mask

    def help(self):
        print(" Add some helps later on ")

    def print_cpu(self, icpu):
        if icpu == max(self.cpus):
            print("Loading particles in {}-th cpu output out of {} cpus.\n"
            .format(icpu, len(self.cpus)))
        else:
            print("Loading particles in {}-th cpu output out of {} cpus.\r"
            .format(icpu, len(self.cpus)))

    def load(self, fortran=True, read_metal=True, verbose=False, ptypes=None, **kwargs):
        """ tests whether the files exist, and then calls load() or load_dmo()
        """
        if verbose:
            self.print_cpu()
        if ptypes is not None: self.set_dtype(ptypes)
        if self.config.sim_type in ["fornax", "nh"]:
            self.load_fortran_new(**kwargs)
        else:
            if self.dmo:
                self.load_dmo(**kwargs)
            else:
                #if fortran:
                return self.load_fortran(read_metal=read_metal, **kwargs)
                #else:
                #    return self.load_general(self, **kwargs)

    def get_dmo_ntot(self):
        ranges = self.ranges
        ndm_tot = 0
        for icpu in self.cpus:
            with open(self._fnbase + str(icpu).zfill(5), "rb") as f:
            # skip header
                header_icpu = read_header(f, self._ramses_particle_header)

                npart_icpu = header_icpu['npart']
                # position
                x_temp = read_fortran(f, np.dtype('f8'), npart_icpu)  # row-major
                y_temp = read_fortran(f, np.dtype('f8'), npart_icpu)
                z_temp = read_fortran(f, np.dtype('f8'), npart_icpu)

                range_ok = np.where((ranges[0][0] < x_temp)
                                    & (ranges[0][1] > x_temp)
                                    & (ranges[1][0] < y_temp)
                                    & (ranges[1][1] > y_temp)
                                    & (ranges[2][0] < z_temp)
                                    & (ranges[2][1] > z_temp))
            ndm_tot = ndm_tot + len(range_ok[0])
        return ndm_tot

    def load_dmo(self, zoom=False, verbose=False, ranges=None):
        """
        DMO run output is simpler:
        no time, no metal, no need to calculate number of each type of particles.
        So it should run faster!
        """
        # function argument is evaluated on function defining time,
        # So you can't pass the actual value of self.info instance
        if ranges is None:
            ranges = self.ranges
        # Total particle number from selected cpus.
        npart_arr = self._get_npart_arr(self.cpus)
        if verbose:
            print("Loading particle... \n ranges:", ranges)
            print(self.cpus)
            print('npart_arr:', npart_arr)
        ndm_tot = self.get_dmo_ntot()

        dtype = self._get_dtype("dm")
        self.dm = np.recarray(ndm_tot, dtype=dtype)
        i_skip_dm = 0

        for icpu in self.cpus:
            if verbose:
                self.print_cpu(icpu)

            with open(self._fnbase + str(icpu).zfill(5), "rb") as f: # +1

                # skip header
                header_icpu = read_header(f, self._ramses_particle_header)

                npart_icpu = header_icpu['npart']

                # position
                x_temp = read_fortran(f, np.dtype('f8'), npart_icpu)  # row-major
                y_temp = read_fortran(f, np.dtype('f8'), npart_icpu)
                z_temp = read_fortran(f, np.dtype('f8'), npart_icpu)

                range_ok = np.where((ranges[0][0] < x_temp)
                                    & (ranges[0][1] > x_temp)
                                    & (ranges[1][0] < y_temp)
                                    & (ranges[1][1] > y_temp)
                                    & (ranges[2][0] < z_temp)
                                    & (ranges[2][1] > z_temp))[0]
                ndm_icpu = len(range_ok)

                self.dm['x'][i_skip_dm:i_skip_dm + ndm_icpu] = x_temp[range_ok]
                self.dm['y'][i_skip_dm:i_skip_dm + ndm_icpu] = y_temp[range_ok]
                self.dm['z'][i_skip_dm:i_skip_dm + ndm_icpu] = z_temp[range_ok]

                # velocity
                if "vx" in self.dm.dtype.names:
                    self.dm['vx'][i_skip_dm:i_skip_dm + ndm_icpu] = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
                    self.dm['vy'][i_skip_dm:i_skip_dm + ndm_icpu] = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
                    self.dm['vz'][i_skip_dm:i_skip_dm + ndm_icpu] = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]

                else:
                    read_fortran(f, np.dtype('f8'), npart_icpu)
                    read_fortran(f, np.dtype('f8'), npart_icpu)
                    read_fortran(f, np.dtype('f8'), npart_icpu)

                # mass
                if "m" in self.dm.dtype.names:
                    m_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
                else:
                    read_fortran(f, np.dtype('f8'), npart_icpu)

                # id
                self.dm['id'][i_skip_dm:i_skip_dm + ndm_icpu] = read_fortran(f, np.dtype('i4'), npart_icpu)[range_ok]

                # refinement
                if self.dm_with_ref:
                    ref_temp = read_fortran(f, np.dtype('i4'), npart_icpu)[range_ok]
                else:
                    read_fortran(f, np.dtype('i4'), npart_icpu)

                if "m" in self.dm.dtype.names:
                    self.dm['m'][i_skip_dm:i_skip_dm + ndm_icpu] = m_temp

                if  "ref" in self.dm.dtype.names:
                    self.dm['ref'][i_skip_dm:i_skip_dm + ndm_icpu] = ref_temp[range_ok]
                i_skip_dm += ndm_icpu

    def _get_dtype(self, ptype):
        part_now = getattr(self.ptypes, ptype)
        dtype=[]
        if part_now.pos:
            dtype.append(('x', '<f8'))
            dtype.append(('y', '<f8'))
            dtype.append(('z', '<f8'))

        if part_now.vel:
            dtype.append(('vx', '<f8'))
            dtype.append(('vy', '<f8'))
            dtype.append(('vz', '<f8'))

        if part_now.mass:
            dtype.append(('m', '<f8'))

        if part_now.id:
            dtype.append(('id', '<i4'))

        if part_now.time:
            dtype.append(('time', '<f8'))

        if part_now.metal:
            dtype.append(('metal', '<f8'))

        return dtype

    def load_fortran_new(self, target_fields=None, cpulist=None):
        """Reads particle data from current box.

        Parameters
        ----------
        target_fields: list of string
            target field name to read. If None is passed, read all fields.
        cpulist: list of integer
            target list of cpus, if specified, read selected cpus regardless of current box.
        Returns
        -------
        part : RamsesSnapshot.Particle object
            particle data, can be accessed as attributes also.

        example
        -------
        ss = pyram.load.sim.Sim(193, base="./", sim_type="fornax")
        ss.add_part(ptypes = ["tracer", "dm", "star", "star_tracer", "sink"])

        """
        from .dtypes import fornax_dtypes
        #from ..config import part_family
        # Move to config. -- There's already all the dicts needed. 
        
        #family_keys = [0,1,2,-2,3]

        #(xmi, xma), (ymi, yma), (zmi, zma) = self.ranges

        if(cpulist is None):
            cpulist = self.info.cpus
        
        if (cpulist.size > 0):
            wdir = self.base + '/snapshots/'

            readr.read_part(wdir, self.nout, cpulist, self.config.sim_type, 
                            np.array(self.info.ranges).ravel(), True, self.config.longint)
            
            if self.config.sim_type == 'nh':
                part_float = fromarrays([*readr.real_table.T],
                 dtype=[ii for ii in self.rur_dtype if ii[1] == 'f8'])
                part_int = fromarrays([*readr.integer_table.T],
                 dtype=[('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')])
            else:
                #timer.record()
                # If outsdie ROI, byte_table(:,1)=-128
                if(self.config.longint):
                    arr = [*readr.real_table.T, 
                        *readr.long_table.T 
                        *readr.integer_table.T, 
                        *readr.byte_table.T]
                else:
                    arr = [*readr.real_table.T, 
                        *readr.integer_table.T, 
                        *readr.byte_table.T]

            #timer.start('Building table for %d particles... ' % readr.integer_table.shape[0], 1)
            #part = fromarrays(arr, dtype=fornax_dtypes['raw_dtype'])[ind_ok]
            if self.config.sim_type == 'nh':
                from .dtypes import nh_dtypes
                #print(part_float.dtype)
                if self.nstar > 0:
                    isnt_star = part_float['time'] == 0
                    
                    if 'star' in self.pt:
                        dtype_star = nh_dtypes['dtype_star']

                        #if read_metal:
                        dtype_star.update({'metal': (('<f8', 1), 72)})
                        istar = np.where(~isnt_star)[0]
                        self.star = np.zeros(len(istar), dtype=dtype_star)
                        for i, tag in enumerate(['x','y','z','vx','vy','vz','m', 'time', 'metal']):
                            self.star[tag] = part_float[tag][istar]
                        #if read_metal:
                        #self.star['metal'] = part_float[tag][istar]
                        self.star['id'] = part_int['id'][istar]
                    if 'sink' in self.pt:
                        dtype_sink = dtype_dm
                        is_psuedo = isnt_star * (part_int < 0)
                        ind_sink = np.where(is_psuedo*(part_float[:,6] > 0))[0]
                        self.sink = np.zeros(self.nsink, dtype=dtype_sink)
                        for i, tag in enumerate(['x','y','z','vx','vy','vz','m']):
                            self.sink[tag] = part_float[ind_sink,i]
                        self.sink['id'] = part_int['id'][ind_sink]
                else:
                    isnt_star = np.arange(len(part_float))
                if 'dm' in self.pt:
                    dtype_dm = nh_dtypes['dtype_dm']
                    ind_dm = np.where(isnt_star * part_int['id'] > 0)[0]
                    self.dm = np.zeros(len(ind_dm), dtype=dtype_dm)
                    for i, tag in enumerate(['x','y','z','vx','vy','vz','m']):
                        self.dm[tag] = part_float[tag][ind_dm]
                    self.dm['id'] = part_int['id'][ind_dm]
                if 'tracer' in self.pt:
                    """
                    Todo 

                    implement tracer reader - or is it in part?
                    """
                    dtype_tracer = nh_dtypes["dtype_tracer"]
                    ind_tracer = np.where(is_psuedo * (part_float[:,6] == 0))[0]
                    self.tracer = np.zeros(self.ntracer, dtype=dtype_tracer)
                    for i, tag in enumerate(['x','y','z','vx','vy','vz']):
                        self.tracer[tag] = part_float[ind_tracer,i]
                    self.tracer['id'] = part_int[ind_tracer]
            else:
                bt = fromarrays([*readr.byte_table.T], dtype=[("family",'i1'),('tag','i1')])
                ind_ok = np.where(bt['family'] > -128)[0]
                part = fromarrays(arr, dtype=self.rur_dtype)[ind_ok]
                
                # copy memory
                for (pt, fam) in zip(self.pt, self.family_keys):
                    #print('family', fam)
                    tmp = part[part['family']==fam]
                    #print("FORNAX dtype", fornax_dtypes[pt])
                    to_store = np.zeros(len(tmp), dtype=fornax_dtypes[pt])
                    for dt in fornax_dtypes[pt]:
                        try:
                            to_store[dt] = tmp[dt]
                        except:
                            print("skipping", dt)
                    #print("setting attributes:", pt)
                    setattr(self, pt, to_store)
            
            # deallocate fortran memory
            readr.clear_all()
            

    def load_fortran(self,
                     return_meta=False,
                     read_metal=True,
                     return_whole=False,
                     verbose=False):
        """
        # parameters:

        # NOTE:
        part_load_module.load_part returns only two types of particles: star and DM.
        It's Python's job to separate particles into DM, sink, and tracer.

        # distinguish sink / dm / star
        # non star : t == 0
        # DM : t == 0, id > 0
        # sink : t ==0, id < 0, mass > 0
        # tracer : t ==0, id < 0, mass = 0

        Note on fortran array allocation & deallocation.
        
        < m.f90 >
        module m
          real, allocatable, dimension(:) :: x
          contains
            subroutine foo
              ...
          ...
        <In Python session>
        >>> m.foo() # call F90 module m function foo
        >>> m.x = [1, 2, 3] # allocate F90 x and initialize with given values
        >>> m.x # returns the content of x as array
        >>> m.x = None # deallocate F90 x array

        Or, Like what San did -- readr.clear_all() --,
        make a Fortran subroutine to deallocate all allocated arrays in Fortran.

        """
        if verbose: print("Loading by fortran module")
        
        (xmi, xma), (ymi, yma), (zmi, zma) = self.ranges
        work_dir = self.base + '/snapshots/output_' + str(self.nout).zfill(5)
        self.ndm, self.nstar, self.nsink, self.ntracer = part_load_module.count_part( \
                            work_dir, xmi, xma, ymi, yma, zmi, zma, self.cpus)
        npart_tot = np.int32(sum((self.ndm, self.nstar, self.nsink, self.ntracer)))
        if verbose:
            print("star:", self.nstar, "DM:", self.ndm,
                  "sink:", self.nsink, "tracer:",self.ntracer)
        if npart_tot == 0:
            return

        # I want to pass shared arrays that are allocated in Python side.
        # But I get the following message :(
        # ValueError: failed to initialize intent(inout) array -- input not fortran contiguous
        # ndm_actual = max([self.ndm, 1]) # why 1?
        # nstar_actual = max([self.nstar, 1])
        read_metal = 1 #
        part_float = np.zeros((npart_tot,9), order="F",dtype=np.float64)
        part_int = np.zeros(npart_tot, order="F", dtype=np.int32)
        part_load_module.load_part(part_float, part_int,
                npart_tot,
                work_dir, xmi, xma, ymi, yma, zmi, zma, read_metal, self.cpus)

        isnt_star = part_float[:,7] == 0
        from .dtypes import nh_dtypes
        if 'star' in self.pt:
            dtype_star = nh_dtypes['dtype_star']

            if read_metal:
                dtype_star.update({'metal': (('<f8', 1), 72)})
            istar = np.where(~isnt_star)[0]
            self.star = np.zeros(self.nstar, dtype=dtype_star)
            for i, tag in enumerate(['x','y','z','vx','vy','vz','m', 'time']):
                self.star[tag] = part_float[istar,i]
            if read_metal:
                self.star['metal'] = part_float[istar,8]
            self.star['id'] = part_int[istar]
        if 'dm' in self.pt:
            dtype_dm = nh_dtypes['dtype_dm']
            ind_dm = np.where(isnt_star * part_int > 0)[0]
            self.dm = np.zeros(self.ndm, dtype=dtype_dm)
            for i, tag in enumerate(['x','y','z','vx','vy','vz','m']):
                self.dm[tag] = part_float[ind_dm,i]
            self.dm['id'] = part_int[ind_dm]
        if 'sink' in self.pt:
            dtype_sink = dtype_dm
            is_psuedo = isnt_star * (part_int < 0)
            ind_sink = np.where(is_psuedo*(part_float[:,6] > 0))[0]
            self.sink = np.zeros(self.nsink, dtype=dtype_sink)
            for i, tag in enumerate(['x','y','z','vx','vy','vz','m']):
                self.sink[tag] = part_float[ind_sink,i]
            self.sink['id'] = part_int[ind_sink]
        if 'tracer' in self.pt:
            """
            Todo 

            implement tracer reader - or is it in part?
            """
            dtype_tracer = nh_dtypes["dtype_tracer"]
            ind_tracer = np.where(is_psuedo * (part_float[:,6] == 0))[0]
            self.tracer = np.zeros(self.ntracer, dtype=dtype_tracer)
            for i, tag in enumerate(['x','y','z','vx','vy','vz']):
                self.tracer[tag] = part_float[ind_tracer,i]
            self.tracer['id'] = part_int[ind_tracer]

        if verbose: print("Fortran-reading done")

        part_float=0
        part_int=0
        if return_whole:
            return part_float, part_int
        # now, how to deallocate it?

    def reload(self, ranges=None, verbose=False):
        """
        Returns sub array of current particle data.
        """
        pass
        # fancy indexing returns view to array.
        # BUT indexing to multiple fields does not return a view, but copies memory.

    def load_sink(self, icoarse=None, read_props=True, read_raw=True):
        self._sink = Sink(self.info, self.config)
        if read_props:
            self.sinkp = self._sink.read_sinkprop( self._sink._path, icoarse=icoarse)
        if read_raw:
            self.sinkr = self._sink.read_sink_raw(1)

class Sink(Simbase):
    def __init__(self, info, config, path='./'):
        self.info = info
        self.config = config
        try:
            path = self.config.sink_path
        except:
            pass
        self._path = path

    def read_sinkprop(self, 
                    path_in_repo='', 
                    icoarse=None, 
                    drag_part=True, 
                    raw_data=False, 
                    return_aexp=False
                    ):
        if(icoarse is None):
            icoarse = self.info.nstep_coarse-1
        
        path = join(self.info.base, path_in_repo)
        check = join(path, sinkprop_format.format(icoarse=icoarse))
        if(not exists(check)):
            raise FileNotFoundError('Sinkprop file not found: %s' % check)
        

        readr.read_sinkprop(path, icoarse, drag_part, self.config['sim_type'])
        arr = [*readr.integer_table.T, *readr.real_table.T]

        #timer.start('Building table for %d smbhs... ' % readr.integer_table.shape[0], 1)
        if(raw_data):
            return arr
        if(drag_part):
            dtype = sink_prop_dtype_drag
        else:
            dtype = sink_prop_dtype
        if(self.config.sim_type == 'fornax'):
            dtype = sink_prop_dtype_drag_fornax
        if(len(arr) != len(dtype)):
            readr.clear_all()
            raise ValueError('Number of fields mismatch\n'
                            'Recieved: %d, Allocated: %d' % (len(arr), len(dtype)))
        sink = fromarrays(arr, dtype=dtype)
        #timer.record()
        aexp = np.copy(readr.aexp)
        readr.clear_all()

        if(return_aexp):
            return sink, aexp
        else:
            return sink

    def read_sink_raw(self, icpu, path_in_repo=''):
        path = join(self.info.base, path_in_repo)
        readr.read_sink(path, self.info.nout, icpu, self.info.lmin, self.info.lmax)
        arr = [*readr.integer_table.T, *readr.real_table.T]
        sink = fromarrays(arr)
        readr.clear_all()
        return sink
