# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 22:46:56 2015

@author: hoseung
"""


def count_part(part_obj, zoom=False, verbose=False, ranges=None, pq=None):
    import numpy as np
    from load.utils import read_header, read_fortran    
    # only xyz coordinate is useful with Hilbert space domain decomposition
    # information.
    '''
    takes meta data from Part object.
    But, data is stored to the out_array rather than
    being added as a new attribute to the object.
    
    '''
    if ranges is None:
        ranges = part_obj.ranges
    print("Loading particle... \n ranges:", ranges)
    # Total particle number from selected cpus.
    npart_arr = part_obj._get_npart_arr(part_obj.cpus)
#        print(part_obj.info.cpus)
    print('npart_arr:', npart_arr)

    # Calculate total number of DM, star, sink from selected cpus.
    # partilce ID, time are needed to distinguish particle type.
    ndm_tot = 0
    nstar_tot = 0
    nsink_tot = 0

    for icpu in part_obj.cpus:
        if verbose:
            part_obj.print_cpu(icpu)
        with open(part_obj._fbase + str(icpu).zfill(5), "rb") as f: # +1
            header_icpu = read_header(f, part_obj._ramses_particle_header)
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
            # star particles have positive creation time.
            # mostly positive ID, but young stars
            # (before SN burst) have negative IDs.
            i_star = (abs(t_temp) > 0.0000001)
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

#        part_obj.ndm = sum(npart_arr) - part_obj.nstar - part_obj.nsink * 2109

    # Total number of particles stored in the cpus.
    # But particles within ranges, not in cpus, are eventually returned.
    # So ndm_tot is not going to be the size of DM array.
    if 'dm' in part_obj.pt:
#        dtype = part_obj._get_dtype("dm")
        part_obj.ndm=ndm_tot
#        part_obj.dm = np.recarray(ndm_tot, dtype=dtype)

    if 'star' in part_obj.pt:
#        dtype = part_obj._get_dtype("star")
        part_obj.nstar=nstar_tot
#        part_obj.star = np.recarray(nstar_tot, dtype=dtype)

    if 'sink' in part_obj.pt:        
#        dtype = part_obj._get_dtype("sink")
        part_obj.nsink=nsink_tot
#        part_obj.sink = np.recarray(nsink_tot * 2109, dtype=dtype)

    # part_obj.ndm = ndm_tot
    # part_obj.nstar = nstar_tot
    # part_obj.nsink = nsink_tot / 2109

    print("Total DM particle %d" % ndm_tot)
    print("Total star particle %d" % nstar_tot)
    print("Total sink particle %d" % nsink_tot)
    
    return


#%%
def part2arr(part_obj, dm=None, star=None, sink=None,
             verbose=False):

    import numpy as np
    from load.utils import read_header, read_fortran      
    ranges = part_obj.ranges
    print("Loading particle... \n ranges:", ranges)

    if 'dm' in part_obj.pt:
        assert dm.shape[1] == part_obj.ndm, ("given array size is {}"
        ", whereas obj.ndm is {}".format(dm.shape[1], part_obj.ndm))
        # ( str1 str2 str3.format(vars)) == string concatenation
#        dtype = part_obj._get_dtype("dm")
        i_skip_dm = 0

    if 'star' in part_obj.pt:
        assert star.shape[1] == part_obj.nstar, ("given array size is {}" 
        ", whereas obj.nstar is {}".format(star.shape[1], part_obj.nstar))
        i_skip_star = 0

    # Or, hasattr(part_obj, 'sink')
    if 'sink' in part_obj.pt:        
        assert sink.shape[1] == part_obj.nsink, ("given array size is {}"
        ", whereas obj.ndm is {}".format(sink.shape[1], part_obj.nsink))
        i_skip_sink = 0

    for icpu in part_obj.cpus:
        if verbose:
            part_obj.print_cpu(icpu)

        f = open(part_obj._fbase + str(icpu).zfill(5), "rb")  # +1

        header_icpu = read_header(f, part_obj._ramses_particle_header)
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
        if "vel" in part_obj.pqset:
            vx_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
            vy_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
            vz_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
        else:
            for i in range(3):
                read_fortran(f, np.dtype('f8'), npart_icpu)

        # mass
        if "mass" in part_obj.pqset:
            m_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
        else:
            read_fortran(f, np.dtype('f8'), npart_icpu)

        # id
        id_temp = read_fortran(f, np.dtype('i4'), npart_icpu)[range_ok]

        # refinement
        if "ref" in part_obj.pqset:
            ref_temp = read_fortran(f, np.dtype('i4'), npart_icpu)[range_ok]
        else:
            read_fortran(f, np.dtype('i4'), npart_icpu)
        
        # time - necessary 
        t_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]

        # metal
        if "metal" in part_obj.pqset:
            z_temp = read_fortran(f, np.dtype('f8'), npart_icpu)[range_ok]
        else:
            read_fortran(f, np.dtype('f8'), npart_icpu)

        # distinguish sink / dm / star
        # non star : t == 0
        # sink : t ==0, id < 0

# Copy data to form contiguous arrays of particles.
        if 'star' in part_obj.pt:
            i_star = (abs(t_temp) > 0.0000001)
            nstar_icpu = sum(i_star)
            if part_obj.ptypes.star.pos:
                star[0][i_skip_star:i_skip_star + nstar_icpu] = px_temp[i_star]
                star[1][i_skip_star:i_skip_star + nstar_icpu] = py_temp[i_star]
                star[2][i_skip_star:i_skip_star + nstar_icpu] = pz_temp[i_star]
            if part_obj.ptypes.star.vel:
                star[3][i_skip_star:i_skip_star + nstar_icpu] = vx_temp[i_star]
                star[4][i_skip_star:i_skip_star + nstar_icpu] = vy_temp[i_star]
                star[5][i_skip_star:i_skip_star + nstar_icpu] = vz_temp[i_star]
            if part_obj.ptypes.star.mass:
                star[6][i_skip_star:i_skip_star + nstar_icpu] = m_temp[i_star]
            if part_obj.ptypes.star.id:
                star[7][i_skip_star:i_skip_star + nstar_icpu] = id_temp[i_star]
            if part_obj.ptypes.star.time:
                star[8][i_skip_star:i_skip_star + nstar_icpu] = t_temp[i_star]
            if part_obj.ptypes.star.metal:
                star[9][i_skip_star:i_skip_star + nstar_icpu] = z_temp[i_star]
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

        if 'dm' in part_obj.pt:
            if part_obj.ptypes.star.pos:
                dm[0][i_skip_dm:i_skip_dm + ndm_icpu] = px_temp[i_dm]
                dm[1][i_skip_dm:i_skip_dm + ndm_icpu] = py_temp[i_dm]
                dm[2][i_skip_dm:i_skip_dm + ndm_icpu] = pz_temp[i_dm]
            if part_obj.ptypes.dm.vel:
                dm[3][i_skip_dm:i_skip_dm + ndm_icpu] = vx_temp[i_dm]
                dm[4][i_skip_dm:i_skip_dm + ndm_icpu] = vy_temp[i_dm]
                dm[5][i_skip_dm:i_skip_dm + ndm_icpu] = vz_temp[i_dm]
            if part_obj.ptypes.dm.mass:
                dm[6][i_skip_dm:i_skip_dm + ndm_icpu] = m_temp[i_dm]
            if part_obj.ptypes.dm.id:
                dm[7][i_skip_dm:i_skip_dm + ndm_icpu] = id_temp[i_dm]
            i_skip_dm += ndm_icpu

        # Which is faster?
        # i_star[i_dm] as ndm array
        # or i_dm as npart array + i_sink as npart array
        if 'sink' in part_obj.pt:
            if part_obj.ptypes.dm.pos:
                part_obj.sink['x'][i_skip_sink:i_skip_sink + nsink_icpu] = px_temp[i_sink]
                part_obj.sink['y'][i_skip_sink:i_skip_sink + nsink_icpu] = py_temp[i_sink]
                part_obj.sink['z'][i_skip_sink:i_skip_sink + nsink_icpu] = pz_temp[i_sink]
            if part_obj.ptypes.dm.vel:
                part_obj.sink['vx'][i_skip_sink:i_skip_sink + nsink_icpu] = vx_temp[i_sink]
                part_obj.sink['vy'][i_skip_sink:i_skip_sink + nsink_icpu] = vy_temp[i_sink]
                part_obj.sink['vz'][i_skip_sink:i_skip_sink + nsink_icpu] = vz_temp[i_sink]
            if part_obj.ptypes.dm.mass:
                part_obj.sink['m'][i_skip_sink:i_skip_sink + nsink_icpu] = m_temp[i_sink]
            if part_obj.ptypes.dm.id:
                part_obj.sink['id'][i_skip_sink:i_skip_sink + nsink_icpu] = id_temp[i_sink]
            i_skip_sink += nsink_icpu


#%%

import load
import tree
import numpy as np
import utils.sampling as smp
from utils import util

#ncore = int(input("How man cores? \n"))
#wdir = input("Working directory \n")
wdir = '/home/hoseung/Work/data/AGN2/'
ncore = 2
nout=132
snout = str(nout)
rscale = 0.8
npix=800

info = load.info.Info(nout=nout, base=wdir)

frefine= 'refine_params.txt'
fnml = 'cosmo_200.nml'

ptypes=["star id pos mass vel", "dm id pos mass vel"]

# Load all halo
hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="RS", info=info)
#hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="HM", info=info)

hall.load()
# convert to code unit. - done by default
#hall.normalize()

# subset of halos ONLY inside zoom-in region
i_center = np.where(hall.data['np'] == max(hall.data['np']))[0]
h_ind = smp.extract_halos_within(hall, i_center, scale=2.0)
h = tree.halomodule.Halo()
h.derive_from(hall, h_ind)
#h.derive_from(hall, h_ind)#,  [4921, 5281, 5343, 5365, 5375, 5412, 5415], 5442, 5639, 5665, 6095])

region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)
print(region)

hind = np.where((hall.data.x > region["xr"][0]) & (hall.data.x < region["xr"][1]) &
                (hall.data.y > region["yr"][0]) & (hall.data.y < region["yr"][1]) &
                (hall.data.z > region["zr"][0]) & (hall.data.z < region["zr"][1]) &
                (hall.data.mvir > 1e11))[0]
h.derive_from(hall, hind[5:8])
region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)

s = load.sim.Sim()
s.setup(nout, wdir)

s.set_ranges(region["ranges"])
s.show_cpus()
s.add_part(ptypes)

#%%
count_part(s.part, zoom=False, verbose=False, ranges=None, pq=None)

#%%
def static_array2d(shape):
    import multiprocessing
    import ctypes
    import numpy as np
    shared_array_base = multiprocessing.Array(ctypes.c_double, shape[0]*shape[1])
    # Some handles to the memory? 
    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
    # reshape => add convenience index
    return shared_array.reshape(shape[0], shape[1])
#%%
star = static_array2d((10, int(s.part.nstar)))
dm = static_array2d((8, int(s.part.ndm)))
#%%
part2arr(s.part, dm=dm, star=star, sink=None, verbose=True)


