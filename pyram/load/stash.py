"""
Here I keep unused codes, hoping they will be useful at some point in the near future.

"""
# class Part():
    def load_general(self, zoom=False, verbose=False, ranges=None, pq=None):
        # only xyz coordinate is useful with Hilbert space domain decomposition
        # information.
        if ranges is None:
            ranges = self.ranges
        print("Loading particle... \n ranges:", ranges)
        # Total particle number from selected cpus.
        # npart_tot
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

        # self.ndm = sum(npart_arr) - self.nstar - self.nsink * 2109

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
            if self.config['cosmo']:
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

