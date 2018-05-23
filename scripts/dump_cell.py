
# coding: utf-8

# In[1]:
import numpy as np
import time

def radial_cut(xc,yc,zc,r, x,y,z):
    return np.where(np.square(xc-x) +
                    np.square(yc-y) +
                    np.square(zc-z) < r**2)[0]

def dump_cell(cell, fname, nout, gid):
    from utils import io
    import struct
    
    with open(fname, 'wb') as f:
        # Header 
        # nout
        # galaxy id
        # number of cell
        f.write(struct.pack('i', 4))
        f.write(struct.pack('i', nout))
        f.write(struct.pack('i', 4))
        f.write(struct.pack('i', 4))
        f.write(struct.pack('i', gid))
        f.write(struct.pack('i', 4))
        f.write(struct.pack('i', 4))
        f.write(struct.pack('i', len(cell)))
        f.write(struct.pack('i', 4))
        for field in cell.dtype.names:
            io.write_fortran(f, cell[field])

            
def extract_gas(cell_all, gal, s, rscale=4.0):
    import numpy as np
    import load
    region = [[gal['x'] - rscale * gal['r'], gal['x'] + rscale * gal['r']],
              [gal['y'] - rscale * gal['r'], gal['y'] + rscale * gal['r']],
              [gal['z'] - rscale * gal['r'], gal['z'] + rscale * gal['r']]]
    cpus = np.array(s._hilbert_cpulist(s.info, region))
    return _extract_cell(cell_all, cpus,
           gal['x'], gal['y'], gal['z'], rscale * gal['r'])


def _extract_cell(cell_all, cpus, xc, yc, zc, rr,
                  min_gas_density=0, unit="code"):

    cell_cube = cell_all[np.in1d(cell_all['cpu'], cpus)]

    if len(cell_cube) < 500:
        return None
    ind_c = np.where(np.square(cell_cube['x'] - xc) + 
                     np.square(cell_cube['y'] - yc) +
                     np.square(cell_cube['z'] - zc) < rr**2)[0]

    if min_gas_density > 0:
        ind_c = cell_cube['var0'][ind_c] > min_gas_density

    cell = cell_cube[ind_c]#.copy()
    # Match units to the GM galaxy output. 
    # position in Mpc
    # velocity in kms
    # mass in Msun (gal output originaly in 1e11Msun?)
    # Should I convert dx to physical unit too? 
    if unit == "physical":
        cell['x'] = (cell['x'] -0.5) * info.pboxsize
        cell['y'] = (cell['y'] -0.5) * info.pboxsize
        cell['z'] = (cell['z'] -0.5) * info.pboxsize
        cell['dx'] = cell['dx'] * info.pboxsize

    return cell


from analysis.cal_lambda import *  
def main(nout_ini, nout_fi, wdir='./', rcluster_scale = 2.9):
    import load
    import tree.halomodule as hmo
    import numpy as np
    import os
    import pickle
    import analysis.misc
    import utils.match as mtc

    nouts = range(nout_fi, nout_ini -1, -1)

    for nout in nouts:
        print("Nout =", nout)
        s = load.sim.Sim(nout=nout, base=wdir, setup=True)
        s.add_hydro(load=True)
    
        info = load.info.Info(base=wdir, nout=nout)
        gcat = hmo.Halo(base=wdir, is_gal=True, verbose=True, nout=nout)
        # result = analysis.misc.load_cat(wdir + 'easy_new/catalog' + str(nout) + '.pickle')
        result = 
        
        out_dir = wdir + 'GalaxyMaker/CELL_' + str(nout).zfill(5) + '/'
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
    
        i_ok_gals = mtc.match_list_ind(gcat.data['id'], result['id'])

        # small galaxies have small data.
        # Just save them all. 
#        mstar_min = 2 * get_mstar_min(info.aexp)
#        g_ok, halos = get_sample_gal(wdir, nout, info, prg_only_tree, mstar_min)

#        for gal in g_ok.data:
        g_ok = gcat.data[i_ok_gals]
      
        max_dx = np.unique(s.hydro.cell['dx'])[-3]
        # Less than -5. Maybe -4 is OK? but to be conservative..
        all_cell = s.hydro.cell[s.hydro.cell['dx'] <= max_dx]
        for i, gal in enumerate(g_ok):
            print(" {}-th, {}".format(i,gal['id']))
            t0 = time.time()
            gcell = extract_gas(all_cell, gal, s)
            if gcell is None:
                continue
            t1 = time.time()
            print("length :", len(gcell['x']))
            print("extract_gas took", t1 - t0)
            dump_cell(gcell, out_dir + "gal_cells_" + str(gal['id']).zfill(7), nout, gal['id'])
            #print("dump_cell took", time.time() - t1)


if __name__ == '__main__':
    main(88, 173)
