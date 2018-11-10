# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 22:15:03 2015

@author: hoseung
"""
import numpy as np
import load
from utils.io_module import read_header, read_fortran, skip_fortran

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
#
# Part inherits from Simbase, but if I define __init__ once again here,
# the __init__ methods of Simbase and Part does not merge.
# (Think about it, that's very strange).
# Instead, Simbase.__init__ is overridden.
#
#
    def __init__(self, parent=None, nout=None, info=None, dmo=False, ptypes=None, base='./',
                 region=None, ranges=[[0,1]]*3, cpus=None,
                 cpu_fixed=False,
                 data_dir='snapshots/',
                 dmref=False,
                 dmvel=False,
                 dmmass=True,
                 load=False,
                 cosmo=True,
                 fortran=True,
                 verbose=False):
        """
        parameters
        ----------
        parent : a super class.
            If given, all the values of (non callable) attributes from the parent are inherited.
        ptypes : list of particle type and information.
                ["dm id pos"] or ["dm id pos", "star mass vel"]
        dmo : logical
            If True, a faster, DMO read routine invoked (NOT distingushing particle types).
        region : region dict
            only part of snapshot is loaded.
        dmref : logical
            Set True if the snapshot has DM ref information.
        dmvel : logical

        Notes
        -----
        info is required for domain decomposition.

        """
        super(Part, self).__init__()
        self.cosmo = cosmo
        if info is None:
            assert nout is not None, "either info or nout is required"
            from load.info import Info
            info = Info(base=base,nout=nout, cosmo=self.cosmo)
        self.info = info
        self.nout = info.nout
        #self.ptypes = ptypes
        self.cpus = cpus
        self.cpu_fixed=cpu_fixed
        try:
            self.ncpu = len(self.cpus)
        except:
            self.ncpu = 0
        self.nstar = 0
        self.nsink = 0

        if verbose: print("Part", base)
        try:
            self.set_base(info.base)
            if verbose: print(info.base)
        except:
            self.set_base(base)
        self.data_dir = data_dir
        if ptypes is not None: self.setwhattoread(ptypes)
        self.dmo = dmo
        self.dm_with_ref = dmref
        self.dm_with_vel = dmvel
        self.dm_with_mass = dmmass

        self.set_fbase(self.base, data_dir)

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

        if ptypes is not None and load:
            self.load(fortran=fortran)

    def mass2msun(self):
        """
        No use case?
        """
        for ptype in self.pt:
            part = getattr(self, ptype)
            if max(part["m"]) < 100:
                part["m"] *= self.info.msun


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
        snout = str(self.nout).zfill(5)
        self._fbase = path.abspath(path.join(self.base, data_dir +'output_' + snout + '/part_' + snout + '.out'))

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
        self.pqset.intersection_update(quantities)
        # any elements not listed in qunatities are removed.


    def setDmQuantities(self, vel=False, mass=True, ref=False):
        self.dm_with_vel = vel
        self.dm_with_ref = ref

    def _get_basic_info(self):
        try:
            f = open(self._fbase + '00001', "rb")
        except:
            from glob import glob
            parts = glob(self._fbase + "*")
            f = open(parts[0], "rb")

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

    def load(self, fortran=True, read_metal=True, **kwargs):
        """ tests whether the files exist, and then calls load() or load_dmo()
        """
        if self.dmo:
            self.load_dmo(self, **kwargs)
        else:
            if fortran:
                return self.load_fortran(self, read_metal=read_metal, **kwargs)
            else:
                #self.load_2017(self, **kwargs)
                return self.load_general(self, **kwargs)

    def get_dmo_ntot(self):
        ranges = self.ranges
        ndm_tot = 0
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
        if ranges is None:
            ranges = self.ranges
        print("Loading particle... \n ranges:", ranges)
        # Total particle number from selected cpus.
#        npart_tot
        npart_arr = self._get_npart_arr(self.cpus)
        print('npart_arr:', npart_arr)

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


                # Dark Matter particles have 0 creation time, positive ID
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
            self.sink = np.recarray(nsink_tot, dtype=dtype)
            i_skip_sink = 0

        self.ndm = ndm_tot
        self.nstar = nstar_tot
        self.nsink = nsink_tot# / 2109

        print("Total DM particle %d" % ndm_tot)
        print("Total star particle %d" % nstar_tot)
        print("Total sink particle %d (/2109)" % nsink_tot)

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
            if self.cosmo:
                if "metal" in self.pqset:
                    z_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
                else:
                    read_fortran(f, np.dtype('f8'), npart_icpu)

            # distinguish sink / dm / star
            # non star : t == 0
            # sink : t ==0, id < 0
            # tracer : t ==0, id < 0, mass = 0

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

        """
        from load.part_load import part_load_module
        if verbose: print("Loading by fortran module")
        xmi = self.ranges[0][0]
        xma = self.ranges[0][1]
        ymi = self.ranges[1][0]
        yma = self.ranges[1][1]
        zmi = self.ranges[2][0]
        zma = self.ranges[2][1]
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
        #ndm_actual = max([self.ndm, 1]) # why 1?
        #nstar_actual = max([self.nstar, 1])
        read_metal = 1 #
        part_float = np.zeros((npart_tot,9), order="F",dtype=np.float64)
        part_int = np.zeros(npart_tot, order="F", dtype=np.int32)
        part_load_module.load_part(part_float, part_int,
                npart_tot,
                work_dir, xmi, xma, ymi, yma, zmi, zma, read_metal, self.cpus)

        isnt_star = part_float[:,7] == 0
        if 'star' in self.pt:
            dtype_star = { 'id': (('<i8', 1), 0),
                          'pos': (('<f8', (3,)), 8),
                            'x': (('<f8', 1), 8),
                            'y': (('<f8', 1), 16),
                            'z': (('<f8', 1), 24),
                          'vel': (('<f8', (3,)), 32),
                           'vx': (('<f8', 1), 32),
                           'vy': (('<f8', 1), 40),
                           'vz': (('<f8', 1), 48),
                            'm': (('<f8', 1), 56),
                         'time': (('<f8', 1), 64)}
            if read_metal:
                dtype_star.update({'metal': (('<f8', 1), 72)})
            istar = np.where(~isnt_star)[0]
            self.star = np.zeros(self.nstar, dtype=dtype_star)
            self.star['x'] = part_float[istar,0]
            self.star['y'] = part_float[istar,1]
            self.star['z'] = part_float[istar,2]
            self.star['vx'] = part_float[istar,3]
            self.star['vy'] = part_float[istar,4]
            self.star['vz'] = part_float[istar,5]
            self.star['m'] = part_float[istar,6]
            self.star['time'] = part_float[istar,7]
            if read_metal:
                self.star['metal'] = part_float[istar,8]
            self.star['id'] = part_int[istar]

        if 'dm' in self.pt:
            dtype_dm = { 'id': (('<i8', 1), 0),
                         'pos': (('<f8', (3,)), 8),
                            'x': (('<f8', 1), 8),
                            'y': (('<f8', 1), 16),
                            'z': (('<f8', 1), 24),
                          'vel': (('<f8', (3,)), 32),
                           'vx': (('<f8', 1), 32),
                           'vy': (('<f8', 1), 40),
                           'vz': (('<f8', 1), 48),
                            'm': (('<f8', 1), 56)}
            ind_dm = np.where(isnt_star * part_int > 0)[0]
            self.dm = np.zeros(self.ndm, dtype=dtype_dm)
            self.dm['x'] = part_float[ind_dm,0]
            self.dm['y'] = part_float[ind_dm,1]
            self.dm['z'] = part_float[ind_dm,2]
            self.dm['vx'] = part_float[ind_dm,3]
            self.dm['vy'] = part_float[ind_dm,4]
            self.dm['vz'] = part_float[ind_dm,5]
            self.dm['m'] = part_float[ind_dm,6]
            self.dm['id'] = part_int[ind_dm]
        if 'sink' in self.pt:
            dtype_sink = dtype_dm
            is_psuedo = isnt_star * (part_int < 0)
            ind_sink = np.where(is_psuedo*(part_float[:,6] > 0))[0]
            self.sink = np.zeros(self.nsink, dtype=dtype_sink)
            self.sink['x'] = part_float[ind_sink,0]
            self.sink['y'] = part_float[ind_sink,1]
            self.sink['z'] = part_float[ind_sink,2]
            self.sink['vx'] = part_float[ind_sink,3]
            self.sink['vy'] = part_float[ind_sink,4]
            self.sink['vz'] = part_float[ind_sink,5]
            self.sink['m'] = part_float[ind_sink,6]
            self.sink['id'] = part_int[ind_sink]
        if 'tracer' in self.pt:
            dtype_dm = { 'id': (('<i8', 1), 0),
                         'pos': (('<f8', (3,)), 8),
                            'x': (('<f8', 1), 8),
                            'y': (('<f8', 1), 16),
                            'z': (('<f8', 1), 24),
                          'vel': (('<f8', (3,)), 32),
                           'vx': (('<f8', 1), 32),
                           'vy': (('<f8', 1), 40),
                           'vz': (('<f8', 1), 48)}
            ind_tracer = np.where(is_psuedo * (part_float[:,6] == 0))[0]
            self.tracer = np.zeros(self.ntracer, dtype=dtype_trancer)
            self.tracer['x'] = part_float[ind_tracer,0]
            self.tracer['y'] = part_float[ind_tracer,1]
            self.tracer['z'] = part_float[ind_tracer,2]
            self.tracer['vx'] = part_float[ind_tracer,3]
            self.tracer['vy'] = part_float[ind_tracer,4]
            self.tracer['vz'] = part_float[ind_tracer,5]
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
