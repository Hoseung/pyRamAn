
# coding: utf-8

# In[1]:
import numpy as np

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

            
def extract_gas(cell_all, gal, rscale=5.0):
    import numpy as np
    return _extract_cell(cell_all, gal['x'], gal['y'], gal['z'], rscale * gal['r'])


def _extract_cell(cell_all, xc, yc, zc, rr,
                  min_gas_density=0, unit="code"):
    ind_c = np.where((np.square(cell_all['x'] - xc) + 
                      np.square(cell_all['y'] - yc) +
                      np.square(cell_all['z'] - zc)) < rr**2)[0]
    
    if min_gas_density > 0:
        ind_c = ind_c * (cell_all['var0'] > min_gas_density)

    cell = cell_all[ind_c].copy()
    # Match units to the GM galaxy output. 
    # position in Mpc
    # velocity in kms
    # mass in Msun (gal output originaly in 1e11Msun?)
    # Should I convert dx in physical unit too? 
    if unit == "physical":
        cell['x'] = (cell['x'] -0.5) * info.pboxsize
        cell['y'] = (cell['y'] -0.5) * info.pboxsize
        cell['z'] = (cell['z'] -0.5) * info.pboxsize
        cell['dx'] = cell['dx'] * info.pboxsize

    return cell


# In[2]:

def main(nout_ini, nout_fi, wdir='./', rcluster_scale = 2.9):
    import load
    import tree.halomodule as hmo
    import numpy as np
    import os
    
    nouts = range(nout_fi, nout_ini, -1)

# In[15]:
    
    for nout in nouts:
        s = load.sim.Sim(nout=nout, base=wdir, setup=True)
        s.add_hydro(load=True)
    
        info = load.info.Info(base=wdir, nout=nout)
        gcat = hmo.Halo(base=wdir, is_gal=True, verbose=True, nout=nout)
        hcat = hmo.Halo(base=wdir, is_gal=False, nout=nout)
    
        cluster = hcat.data[np.argmax(hcat.data['np'])]
        out_dir = wdir + 'GalaxyMaker/CELL_' + str(nout).zfill(5) + '/'
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
    
        i_ok_gals = radial_cut(cluster['x'], cluster['y'], cluster['z'],
                               rcluster_scale * cluster['rvir'],
                               gcat.data['x'], gcat.data['y'], gcat.data['z'])
    
        g_ok = gcat.data[i_ok_gals]
        for gal in g_ok:
            gcell = extract_gas(s.hydro.cell, gal)
            dump_cell(gcell, out_dir + "gal_cells_" + str(gal['id']).zfill(7), nout, gal['id'])



if __name__ == '__main__':
    main(100, 125)
