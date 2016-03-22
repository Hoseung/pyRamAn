# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 22:15:03 2015

@author: hoseung
"""
import numpy as np
import load
from load.utils import read_header, read_fortran

"""
pt = []
pq = []
pqset = set([])
for pp in ptype:
    pp = pp.lower()
    pt.append(pp.split()[0])
    pq.append(pp.split()[1:])
    pqset.update(pp.split()[1:])
quantities = set(["mass", "id", "vel", "ref", "time", "metal"])
pqset.intersection_update(quantities) # any elements not listed in qunatities are removed.
"""  
# Check if there is dm / star / sink in quantities list
# or mass, id, vel, and so on in ptype list.
# Only exact match (lower or upper case) works. No tolerence for errata.

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

class Part(load.sim.Simbase):
    """ 
    Particle data container supports different types of particles,
    and some meta data.
        
    DM, star, sink
    """
# Rule of thumb: initialize everything in __init__
# Otherwise, a method may fail to find an object to act on.
# This is called 'Object consistency'
    def __init__(self, info, dmo=False, ptypes=None, base='./', 
                 data_dir='snapshots', dmref=False, dmvel=False, dmmass=True):
        """        
        parameters
        ----------
        ptypes : list of particle type and information. 
                ["dm id pos"] or ["dm id pos", "star mass vel"]
        """
        self.info = info
        self.ptypes = ptypes
        self.ncpu = 0
        self.nstar = 0
        self.nsink = 0
        self.set_base(base)
        self.data_dir = data_dir
        self.setwhattoread(ptypes)
        self.dmo = dmo
        self.dm_with_ref = dmref
        self.dm_with_vel = dmvel
        self.dm_with_mass = dmmass
        self.set_ranges(ranges=info.ranges)
#        self.cpus = info.cpus
        
        self.set_fbase(self.base, data_dir)
        # header structure
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

    def mass2msun(self):
        for ptype in self.pt:
            part = getattr(self, ptype)
            if max(part["m"]) < 100:
                part["m"] *= self.info.msun

    def set_ranges(self, ranges=None):
        if ranges is None:
            ranges = self.ranges
        nr = np.asarray(ranges)
        if not(nr.shape[0] == 3 and nr.shape[1] == 2):
            # Because actual operation on the given input(ranges)
            # does not take place soon, it's not a good place to use
            # try & except clause. There is nothing to try yet.
            print(' Error!')
            print('Shape of ranges is wrong:', nr.shape)
            print('example : [[0.1,0.3],[0.2,0.4],[0.6,0.8]] \n')
        else:
            self.ranges = ranges
            self.set_cpus(self._hilbert_cpulist(self.info, self.ranges))

    def set_base(self, base):
        """
            Sets Working directory.
        """
        from os import path
        self.base = path.abspath(base)
#        self.show_base()
        
    def set_fbase(self, base, data_dir):
        """
            Sets Working directory.
        """
        from os import path
        snout = str(self.info.nout).zfill(5)
        self._fbase = path.abspath(path.join(self.base, data_dir +'/output_' + snout + '/part_' + snout + '.out'))

    def set_cpus(self, cpus):
        self.cpus = cpus
        try:
            print("Updating info.cpus")
            self.info._set_cpus(self.get_cpus())
        except AttributeError:
            print("No info._set_cpus attribute??")

    def setwhattoread(self, ptypes):
        """
        Sets which types of quantiies to read.
        Because there is no separation b/w particles, you need to read all ID
        even if you want only DM IDs.
        """
        self.pt = []
        self.pq = []
        self.pqset = set([])
        for pp in ptypes:
            pp = pp.lower()
            self.pt.append(pp.split()[0])
            self.pq.append(pp.split()[1:])
            self.pqset.update(pp.split()[1:])

        self.ptypes = Ptypes(self.pt, self.pq)   
# Check if there is dm / star / sink in quantities list
# or mass, id, vel, and so on in ptype list.
# Only exact match (lower or upper case) works. No tolerence for errata.
        quantities = set(["mass", "id", "vel", "ref", "time", "metal"])
        self.pqset.intersection_update(quantities) # any elements not listed in qunatities are removed.

    def setDmQuantities(self, vel=False, mass=True, ref=False):
        self.dm_with_vel = vel
        self.dm_with_ref = ref

    def _get_basic_info(self):
        f = open(self._fbase + '00001', "rb")
        header = read_header(f, self._ramses_particle_header)
        self.ncpu = header['ncpu']
        self.nstar = header['nstar']
        self.nsink = header['nsink']

    def _get_npart_arr(self, cpus):
        npart_arr = []
        for icpu in cpus:
            with open(self._fbase + str(icpu).zfill(5), "rb") as f:
                header = read_header(f, self._ramses_particle_header)
                npart_arr.append(header['npart'])

        return npart_arr  # ,nstar_this, my_mask

    def search_zoomin(self, scale=1.0, load=False):
        from utils import sampling
        if not hasattr(self, 'dm'):
            if load:
                print("Let's load particle data first!")
                self.load()
            else:
                print("Couldn't find dm data. \n \
                    Make sure you load it first, \n \
                    or use load=True option")

        imin = np.where(self.dm['m'] == self.dm['m'].min())
        xr = [self.dm['x'][imin].min(), self.dm['x'][imin].max()]
        yr = [self.dm['y'][imin].min(), self.dm['y'][imin].max()]
        zr = [self.dm['z'][imin].min(), self.dm['z'][imin].max()]

        xc = 0.5 * sum(xr)
        yc = 0.5 * sum(yr)
        zc = 0.5 * sum(zr)

        radius = 0.5 * max([xr[1]-xr[0], yr[1]-yr[0], zr[1]-zr[0]]) * scale
#        print(radius)

        return(sampling.set_region(centers=[xc, yc, zc], radius=radius))

    def help():
        print(" Add some helps later on ")

    def print_cpu(self, icpu):
        if icpu == max(self.cpus):
            print("Loading particles in {}-th cpu output out of {} cpus.\n"
            .format(icpu, len(self.cpus)))
        else:
            print("Loading particles in {}-th cpu output out of {} cpus.\r"
            .format(icpu, len(self.cpus)))

    def load(self, fortran=False, **kwargs):
        """ tests whether the files exist, and then calls load() or load_dmo()
        """
        if self.dmo:
            self.load_dmo(self, **kwargs)
        else:
            if fortran:
                self.load_fortran(self, **kwargs) 
            else:
                self.load_general(self, **kwargs)

    def get_dmo_ntot(self):
#        ranges = self.info.ranges
        ranges = self.ranges
        ndm_tot = 0
#        for icpu in self.info.cpus:
        for icpu in self.cpus:
            with open(self._fbase + str(icpu).zfill(5), "rb") as f:
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
        DMO run output is much simpler:
        no time, no metal, no need to calculate number of each type of particles.
        So it should run much faster!
        """
        # function argument is evaluated on function defining time,
        # So you can't pass the actual value of self.info instance
        if ranges is None:
            ranges = self.info.ranges
        # Total particle number from selected cpus.
#        npart_arr = self._get_npart_arr(self.info.cpus)
        npart_arr = self._get_npart_arr(self.cpus)
        if verbose:
            print("Loading particle... \n ranges:", ranges)
            print(self.cpus)
            print('npart_arr:', npart_arr)
        ndm_tot = self.get_dmo_ntot()

        dtype = self._get_dtype("dm")
        self.dm = np.recarray(ndm_tot, dtype=dtype)
        i_skip_dm = 0

#        for icpu in self.info.cpus:
        for icpu in self.cpus:
            if verbose:
                self.print_cpu(icpu)

            with open(self._fbase + str(icpu).zfill(5), "rb") as f: # +1
    
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

    def load_general(self, zoom=False, verbose=False, ranges=None, pq=None):
        # only xyz coordinate is useful with Hilbert space domain decomposition
        # information.
        '''
        load(self,zoom=False,xr=[0,1],yr=[0,1],zr=[0,1]):
        self.dm['px'] = np.zeros((self.ndim,self.ndm), dtype='f8')
        self.dm.vel = np.zeros((self.ndim,self.ndm), dtype='f8')
        self.dm.m = np.zeros(self.ndm, dtype='f8')
        # Mass may also be omitted !
        self.dm.id = np.zeros(self.ndm, dtype='f8')
        #self.dm.t = np.zeros(self.ndm, dtype='f8')
        #self.dm.z = np.zeros(self.ndm, dtype='f8') no need to store
        self.star.pos = np.zeros((self.ndim,self.nstar), dtype='f8')
        self.star.vel = np.zeros((self.ndim,self.nstar), dtype='f8')
        self.star.m = np.zeros(self.ndm, dtype='f8')
        self.star.id = np.zeros(self.ndm, dtype='f8')
        self.star.t = np.zeros(self.ndm, dtype='f8')
        self.star.z = np.zeros(self.ndm, dtype='f8')                
        '''
        if ranges is None:
#            ranges = self.info.ranges
            ranges = self.ranges
        print("Loading particle... \n ranges:", ranges)
        # Total particle number from selected cpus.
#        npart_tot
#        npart_arr = self._get_npart_arr(self.info.cpus)
        npart_arr = self._get_npart_arr(self.cpus)
#        print(self.info.cpus)
        print('npart_arr:', npart_arr)
        #npart_tot = sum(npart_arr)

        if hasattr(self.ptypes, "dm"):
            i_skip_dm = 0
        if hasattr(self.ptypes, "star"):
            i_skip_star = 0
        if hasattr(self.ptypes, "sink"):
            i_skip_sink = 0

        # Calculate total number of DM, star, sink from selected cpus.
        # partilce ID, time are needed to distinguish particle type.
        ndm_tot = 0
        nstar_tot = 0
        nsink_tot = 0

        for icpu in self.cpus:
            if verbose:
                self.print_cpu(icpu)
            with open(self._fbase + str(icpu).zfill(5), "rb") as f: # +1
                header_icpu = read_header(f, self._ramses_particle_header)
                npart_icpu = header_icpu['npart']
                if verbose:
                	print('cpu %s has %s particles.' % (icpu, npart_icpu))
                
                # read position and determine number of particles in the ranges.
                x_temp = read_fortran(f, np.dtype('f8'), npart_icpu)  # row-major
                y_temp = read_fortran(f, np.dtype('f8'), npart_icpu)
                z_temp = read_fortran(f, np.dtype('f8'), npart_icpu)
                
                range_ok = np.where((ranges[0][0] < x_temp)
                					& (ranges[0][1] > x_temp)
                					& (ranges[1][0] < y_temp)
                					& (ranges[1][1] > y_temp)
                					& (ranges[2][0] < z_temp)
                					& (ranges[2][1] > z_temp))
                
                # skip velocity
                read_fortran(f, np.dtype('f8'), npart_icpu)
                read_fortran(f, np.dtype('f8'), npart_icpu)
                read_fortran(f, np.dtype('f8'), npart_icpu)
                
                # skip mass
                read_fortran(f, np.dtype('f8'), npart_icpu)
                
                # read particle id
                id_temp = read_fortran(f, np.dtype('i4'), npart_icpu)[range_ok]
                
                # skip refinement
                read_fortran(f, np.dtype('i4'), npart_icpu)
                
                # read time
                t_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
                # star particles have positive creation time.. no!
                # mostly positive ID, but young stars
                # (before SN burst) have negative IDs.
                i_star = (abs(t_temp) != 0.0) # negative creation time for young star!
                i_dm = np.logical_and(id_temp > 0, t_temp == 0)
                # Sink (BH) particles have negative ID and creation time 0.
                i_sink = np.logical_and(id_temp < 0, t_temp == 0)
                
                # Dark Matter particles have 0 creatino time, positive ID
                nstar_icpu = sum(i_star)
                ndm_icpu = sum(i_dm)
                nsink_icpu = sum(i_sink)
                ndm_tot += ndm_icpu
                nstar_tot += nstar_icpu
                nsink_tot += nsink_icpu

        # number of darkmatter = npart - nstar - nsink * 2109
        # But!! nstar and nsink is for the whole simulation volume while
        # npart is for only selected cpus. hmm.

#        self.ndm = sum(npart_arr) - self.nstar - self.nsink * 2109

        # Total number of particles stored in the cpus.
        # But particles within ranges, not in cpus, are eventually returned.
        # So ndm_tot is not going to be the size of DM array.
        if 'dm' in self.pt:
            dtype = self._get_dtype("dm")
            self.dm = np.recarray(ndm_tot, dtype=dtype)
            i_skip_dm = 0

        if 'star' in self.pt:
            dtype = self._get_dtype("star")
            self.star = np.recarray(nstar_tot, dtype=dtype)
            i_skip_star = 0

        # Or, hasattr(self, 'sink')
        if 'sink' in self.pt:        
            dtype = self._get_dtype("sink")
            self.sink = np.recarray(nsink_tot * 2109, dtype=dtype)
            i_skip_sink = 0

        self.ndm = ndm_tot
        self.nstar = nstar_tot
        self.nsink = nsink_tot# / 2109

        print("Total DM particle %d" % ndm_tot)
        print("Total star particle %d" % nstar_tot)
        print("Total sink particle %d" % nsink_tot)


        # iterate over files to read in data
        for icpu in self.cpus:
            if verbose:
                self.print_cpu(icpu)

            f = open(self._fbase + str(icpu).zfill(5), "rb")  # +1

            header_icpu = read_header(f, self._ramses_particle_header)
            # skip header

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

            # make views to the original array
            px_temp = x_temp[range_ok]
            py_temp = y_temp[range_ok]
            pz_temp = z_temp[range_ok]

            # velocity
            if "vel" in self.pqset:
                vx_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
                vy_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
                vz_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
            else:
                for i in range(3):
                    read_fortran(f, np.dtype('f8'), npart_icpu)

            # mass
            if "mass" in self.pqset:
                m_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
            else:
                read_fortran(f, np.dtype('f8'), npart_icpu)

            # id
            id_temp = read_fortran(f, np.dtype('i4'), npart_icpu)[range_ok]

            # refinement
            if "ref" in self.pqset:
                ref_temp = read_fortran(f, np.dtype('i4'), npart_icpu)[range_ok]
            else:
                read_fortran(f, np.dtype('i4'), npart_icpu)
            
            # time - necessary 
            t_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]

            # metal
            if "metal" in self.pqset:
                z_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
            else:
                read_fortran(f, np.dtype('f8'), npart_icpu)

            # distinguish sink / dm / star
            # non star : t == 0
            # sink : t ==0, id < 0

# Copy data to form contiguous arrays of particles.
            if 'star' in self.pt:
                i_star = (abs(t_temp) > 0.0000001)
                nstar_icpu = sum(i_star)
                if self.ptypes.star.pos:
                    self.star['x'][i_skip_star:i_skip_star + nstar_icpu] = px_temp[i_star]
                    self.star['y'][i_skip_star:i_skip_star + nstar_icpu] = py_temp[i_star]
                    self.star['z'][i_skip_star:i_skip_star + nstar_icpu] = pz_temp[i_star]
                if self.ptypes.star.vel:
                    self.star['vx'][i_skip_star:i_skip_star + nstar_icpu] = vx_temp[i_star]
                    self.star['vy'][i_skip_star:i_skip_star + nstar_icpu] = vy_temp[i_star]
                    self.star['vz'][i_skip_star:i_skip_star + nstar_icpu] = vz_temp[i_star]
                if self.ptypes.star.mass:
                    self.star['m' ][i_skip_star:i_skip_star + nstar_icpu] = m_temp[i_star]
                if self.ptypes.star.id:
                    self.star['id'][i_skip_star:i_skip_star + nstar_icpu] = id_temp[i_star]
                if self.ptypes.star.time:
                    self.star['time'][i_skip_star:i_skip_star + nstar_icpu] = t_temp[i_star]
                if self.ptypes.star.metal:
                    self.star['metal'][i_skip_star:i_skip_star + nstar_icpu] = z_temp[i_star]
                i_skip_star += nstar_icpu

            # i_dm = id_temp < 0
            i_dm = np.logical_and(id_temp > 0, t_temp == 0)
            i_sink = np.logical_and(id_temp < 0, t_temp == 0)
            ndm_icpu = sum(i_dm)

            nsink_icpu = sum(i_sink)
            # print('nDM, nSink', ndm_icpu, nsink_icpu)

# Note that if it's two-division separation,
# i_dm = t_temp == 0 and then,
# it's faster to use ~i_dm than to generate another index array.

            if 'dm' in self.pt:
                if self.ptypes.dm.pos:
                    self.dm['x'][i_skip_dm:i_skip_dm + ndm_icpu] = px_temp[i_dm]
                    self.dm['y'][i_skip_dm:i_skip_dm + ndm_icpu] = py_temp[i_dm]
                    self.dm['z'][i_skip_dm:i_skip_dm + ndm_icpu] = pz_temp[i_dm]
                if self.ptypes.dm.vel:
                    self.dm['vx'][i_skip_dm:i_skip_dm + ndm_icpu] = vx_temp[i_dm]
                    self.dm['vy'][i_skip_dm:i_skip_dm + ndm_icpu] = vy_temp[i_dm]
                    self.dm['vz'][i_skip_dm:i_skip_dm + ndm_icpu] = vz_temp[i_dm]
                if self.ptypes.dm.mass:
                    self.dm['m'][i_skip_dm:i_skip_dm + ndm_icpu] = m_temp[i_dm]
                if self.ptypes.dm.id:
                    self.dm['id'][i_skip_dm:i_skip_dm + ndm_icpu] = id_temp[i_dm]
                if self.ptypes.dm.ref:
                    self.dm['ref'][i_skip_dm:i_skip_dm + ndm_icpu] = ref_temp[i_dm]
                i_skip_dm += ndm_icpu

            # Which is faster?
            # i_star[i_dm] as ndm array
            # or i_dm as npart array + i_sink as npart array
            if 'sink' in self.pt:
                if self.ptypes.dm.pos:
                    self.sink['x'][i_skip_sink:i_skip_sink + nsink_icpu] = px_temp[i_sink]
                    self.sink['y'][i_skip_sink:i_skip_sink + nsink_icpu] = py_temp[i_sink]
                    self.sink['z'][i_skip_sink:i_skip_sink + nsink_icpu] = pz_temp[i_sink]
                if self.ptypes.dm.vel:
                    self.sink['vx'][i_skip_sink:i_skip_sink + nsink_icpu] = vx_temp[i_sink]
                    self.sink['vy'][i_skip_sink:i_skip_sink + nsink_icpu] = vy_temp[i_sink]
                    self.sink['vz'][i_skip_sink:i_skip_sink + nsink_icpu] = vz_temp[i_sink]
                if self.ptypes.dm.mass:
                    self.sink['m'][i_skip_sink:i_skip_sink + nsink_icpu] = m_temp[i_sink]
                if self.ptypes.dm.id:
                    self.sink['id'][i_skip_sink:i_skip_sink + nsink_icpu] = id_temp[i_sink]
                i_skip_sink += nsink_icpu

    def load_fortran(self, return_meta=False):
        from load import part_shared
        print("Loading by fortran module")
        xmi = self.info.ranges[0][0]
        xma = self.info.ranges[0][1]
        ymi = self.info.ranges[1][0]
        yma = self.info.ranges[1][1]
        zmi = self.info.ranges[2][0]
        zma = self.info.ranges[2][1]
        work_dir = self.info.base + '/snapshots/output_' + str(self.info.nout).zfill(5)
        
        ndm_actual, nstar_actual, nsink_actual = part_shared.count_part( \
                            work_dir, xmi, xma, ymi, yma, zmi, zma)
#        print(ndm_actual, nstar_actual, nsink_actual)
        self.ndm = ndm_actual
        self.nstar = nstar_actual
        self.nsink = nsink_actual
        if self.ndm == 0 or self.nstar == 0:
            return 
        """
        if return_meta is True:
            return (ndm_actual, nstar_actual, nsink_actual, work_dir, xmi, xma, ymi, yma, zmi, zma)
        else:
        """
# I want to pass shared arrays that are allocated in Python side. 
# But I get the following message
# ValueError: failed to initialize intent(inout) array -- input not fortran contiguous
        ndm_actual = max([ndm_actual, 1])
        nstar_actual = max([nstar_actual, 1])
        star_float, star_int, dm_float, dm_int = part_shared.load_part(
                nstar_actual, ndm_actual, nsink_actual,
                work_dir, xmi, xma, ymi, yma, zmi, zma)

        dtype_star = [('x', '<f8'), ('y', '<f8'), ('z', '<f8'),
                      ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),
                        ('m', '<f8'), ('time', '<f8'), ('metal', '<f8'),
                        ('id', '<i4')]

        self.star = np.zeros(self.nstar, dtype=dtype_star)
        self.star['x'] = star_float[:,0]
        self.star['y'] = star_float[:,1]
        self.star['z'] = star_float[:,2]
        self.star['vx'] = star_float[:,3]
        self.star['vy'] = star_float[:,4]
        self.star['vz'] = star_float[:,5]
        self.star['m'] = star_float[:,6]
        self.star['time'] = star_float[:,7]
        self.star['metal'] = star_float[:,8]
        self.star['id'] = star_int[:]

        dtype_dm = [('x', '<f8'), ('y', '<f8'), ('z', '<f8'),('vx', '<f8'),
                    ('vy', '<f8'), ('vz', '<f8'), ('m', '<f8'), ('id', '<i4')]

        self.dm = np.zeros(self.ndm + self.nsink, dtype=dtype_dm)
        self.dm['x'] = dm_float[:,0]
        self.dm['y'] = dm_float[:,1]
        self.dm['z'] = dm_float[:,2]
        self.dm['vx'] = dm_float[:,3]
        self.dm['vy'] = dm_float[:,4]
        self.dm['vz'] = dm_float[:,5]
        self.dm['m'] = dm_float[:,6]
        self.dm['id'] = dm_int[:]

        print("Fortran-reading done")
        # now, how to deallocate it?

    def reload(self, ranges=None, verbose=False):
        """
        Returns sub array of current particle data.
        """
        pass
        # fancy indexing returns view to array.
        # BUT indexing to multiple fields does not return a view, but copies memory.



