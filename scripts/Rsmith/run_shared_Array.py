# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 17:13:14 2015

@author: hoseung
"""
import numpy as np

def extract_halos_within(halos, ind_center, scale=1.0, Mcut=1e5):
    import numpy as np
    import utils.sampling as smp
    '''
    Returns indices of halos within SCALE * Rvir of the central halo.

    def extract_halos_within(halos, ind_center, scale=1.0)
    halos : halo finder output (single snapshot)
    ind_center : index of central halo
    scale : multiplying factor to the Rvir of the central halo
    '''
    xc = halos['x'][i_center]
    yc = halos['y'][i_center]
    zc = halos['z'][i_center]
    rvir= halos['rvir'][i_center]

    xx = halos['x']
    yy = halos['y']
    zz = halos['z']
    m = np.array(halos['m'])

    dd = smp.distance_to([xc,yc,zc],[xx,yy,zz])

    i_m = m > Mcut
    i_ok = np.logical_and(dd < (rvir * scale), i_m)

    return i_ok
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




def get_idx(tree, hnus, nout=None):
    i_nout = np.where(tree.field('NOUT') == nout)
    i_halo = match_list_ind(tree[i_nout].field('HALNUM'), hnus)

    return tree[i_nout[i_halo]].field('IDX')


def filter_halo_mass(data, Mcut=None):
    m = np.array(data['m'][0])
    #ind = m > Mcut
    #print("# of halos:",sum(ind))
    ind =np.where(m > Mcut)[0]
    print("# of halos:",len(ind))
    return ind


def n_most_massive(data, mass_count=1000):
    m = np.array(data['m'])
    i = np.argsort(m)
    ind = i[:-1 - mass_count:-1]
    return ind


def filter_halo_pnum(data, Ncut=1000):
    npart = np.array(data['np'])
    ind =np.where(npart > Ncut)[0]
    print("# of halos:",len(ind))
    return ind


#%%
''' Cluster 05101, cluster subhaloes (at the final snapshot)
'''

import tree.halomodule as hmo 
from utils import match
import load

fixed_position = True
work_dir = './'
nout_ini = 37
# 27 : z=4
# 37 : z=3
# 20 : ~ z=5
nout_fi = 187
rvir=3.0
# nout_halo = 122 == nout 10
# nout_halo = 0   == nout 132
nouts = range(nout_fi, nout_ini -1, -1) 
Nnouts = len(nouts)

from tree import tmtree
tree = tmtree.load(work_dir=work_dir, filename="halo/TMtree.fits")
tfin = tree[np.where(tree['NOUT'] == 0)]
tini = tree[np.where(tree['NOUT'] == nout_fi - nout_ini)]
#%%
info = load.info.Info()
info.setup(nout=nout_fi, base=work_dir)
info.load()
hh = hmo.Halo(base=work_dir, nout=nout_fi, halofinder='HM', info=info)
hh.load()
halo = hh.data
#halo = hmu.load_data(nout_fi, work_dir=work_dir, normalize=True)
i_center = np.where(halo['np'] == max(halo['np']))
i_satellites = extract_halos_within(halo, i_center, scale=3.0)[0]
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))

# halos found inside the cluster and has tree back to nout_ini
halo_list = halo['id'][i_satellites]
#print(halo_list)
h_ind_ok, halo_ok = tmtree.check_tree_complete(tree, 0, nout_fi - nout_ini, halo_list)
print(len(halo_ok))

try:
    f = open(work_dir + 'satellite_halos.txt', 'w')
except:
    print("No filename is given.\n Try write_halo_xyz(x,y,z,r,filename = fn)")

f.write(" #      ID        x          y         z[Mpc]       vx      vy     vz[km/s]")
f.write("    Rvir(Mpc)      Mvir      Mass[10e10Msun]\n")
zred=[]
# It works, but too slow!

import utils.sampling as smp
from load import part_shared
import time
import multiprocessing as mp
import pandas as pd
import pickle
  
rscale = 0.8
npix=800
rscale_lambda = 3.0
npix_lambda = 30  

ptypes=["star id pos mass vel", "dm id pos mass vel"]
ncore = 2

for inout, nout in enumerate(nouts):
    print(nout)
    info = load.info.Info()
    info.setup(nout=nout, base=work_dir)
    info.load()

    fn = work_dir + 'halo/tree_bricks' + str(nout).zfill(3)
    hall = hmo.Halo(base=work_dir, nout=nout, halofinder='HM', info=info)
    hall.load()

    h_ind = match.match_list_ind(hall.data['id'], halo_ok[:,inout])
    h = hmo.Halo()
    h.derive_from(hall, h_ind)
    
    snout = str(nout)
    fcat = work_dir +"catalog" + snout + ".pickle"
#    galaxy_plot_dir = wdir + 'galaxy_plot' + snout + '/'
        
    region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)   
    s = load.sim.Sim()
    s.setup(nout,  work_dir)
    s.set_ranges(region["ranges"])
    s.show_cpus()
    s.add_part(ptypes)
    s.add_hydro()  
    
    # I assume there is no need to free these shared arrays.
    
    # Set ranges
    xmi = s.info.ranges[0][0]
    xma = s.info.ranges[0][1]
    ymi = s.info.ranges[1][0]
    yma = s.info.ranges[1][1]
    zmi = s.info.ranges[2][0]
    zma = s.info.ranges[2][1]
    work_dir_output = s.info.base + '/snapshots/output_' + str(s.info.nout).zfill(5)
    
    ndm_actual, nstar_actual, nsink_actual = part_shared.count_part(work_dir_output, xmi, xma, ymi, yma, zmi, zma)   
    print(nstar_actual, nsink_actual, ndm_actual)
    
    star_tmp, dm_tmp = part_shared.load_part(
                     nstar_actual, ndm_actual, nsink_actual,
                     work_dir_output, xmi, xma, ymi, yma, zmi, zma)
   
    star = static_array2d((10,nstar_actual))
    dm = static_array2d((8,ndm_actual + nsink_actual)) # no distin
    # pos 3, vel 3, mass, id, time, metal
    # convert to km/s    
    star[0] = star_tmp[:,0]
    star[1] = star_tmp[:,1]
    star[2] = star_tmp[:,2]
    star[3] = star_tmp[:,3] * s.info.kms
    star[4] = star_tmp[:,4] * s.info.kms
    star[5] = star_tmp[:,5] * s.info.kms
    star[6] = star_tmp[:,6] * s.info.msun
    star[7] = star_tmp[:,7]
    star[8] = star_tmp[:,8]
    star[9] = star_tmp[:,9]
    
    if 'dm' in s.part.pt:
        dm[0] = dm_tmp[:,0]
        dm[1] = dm_tmp[:,1]
        dm[2] = dm_tmp[:,2]
        dm[3] = dm_tmp[:,3] * s.info.kms
        dm[4] = dm_tmp[:,4] * s.info.kms
        dm[5] = dm_tmp[:,5] * s.info.kms
        dm[6] = dm_tmp[:,6] * s.info.msun
        dm[7] = dm_tmp[:,7]
    
    s.hydro.amr2cell(lmax=19)
    cell = static_array2d((10,int(s.hydro.cell.shape[0])))
    
    cell[0] = s.hydro.cell['x']
    cell[1] = s.hydro.cell['y']
    cell[2] = s.hydro.cell['z']
    cell[3] = s.hydro.cell['dx']
    cell[4] = s.hydro.cell['var0']
    cell[5] = s.hydro.cell['var1']
    cell[6] = s.hydro.cell['var2']
    cell[7] = s.hydro.cell['var3']
    cell[8] = s.hydro.cell['var4']
    cell[9] = s.hydro.cell['var5']

    def mk_gal(halodata, out_q, info, i, wdir,
                save=False, rscale=0.3, verbose=True, def_param=(dm, star, cell)):
        from galaxymodule import galaxy
    
    #    print(halodata['id'])
        print("{}-th **********************************************".format(i))
    
        gal = galaxy.Galaxy(halodata, radius_method='simple', info=info)
    #    print(gal.halo)
        gal.mk_gal_mp(star=star, dm=dm, cell=cell,
                   rscale=rscale, verbose=True)
        gal_out = {"id":0, "xc":0.0, "yc":0.0, "zc":0.0,
               "vx":0.0, "vy":0.0, "vz":0.0,
               "mstar":0.0, "nstar":0.0, "mgas":0.0,
               "lambda_arr":[], "lambda_r":0, "rgal":0}
    
        if gal.star is False:
            print(gal.id, " Not a good galaxy")
            out_q.put(gal_out)
        else:
            print("R eff:", gal.reff * info.pboxsize)
            gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    
    #        gal.plot_gal(save_dir= galaxy_plot_dir, ioff=True)
        #            gal.save_gal(base=wdir)
    
            # Instead of galaxy class, save them in a dict.    
            gal_out['mstar'] = gal.mstar
            gal_out['mgas'] = gal.mgas
            gal_out['nstar'] = gal.nstar
            gal_out['id'] = gal.id
            gal_out['xc'] = gal.xc * info.pboxsize
            gal_out['yc'] = gal.yc * info.pboxsize
            gal_out['zc'] = gal.zc * info.pboxsize
            gal_out['vx'] = gal.vxc
            gal_out['vy'] = gal.vyc
            gal_out['vz'] = gal.vzc        
            gal_out['lambda_arr'] = gal.lambda_arr
            gal_out['lambda_r'] = gal.lambda_r
            gal_out['rgal'] = gal.reff * info.pboxsize * 1000.0 # in kpc
            out_q.put(gal_out)
        print("{}-th ********************** DONE *************".format(i))
        return True # A return value is expected by Pool.


    # Why cells are in physical unit from the debining whereas part is in code unit? 
    t1 = time.time()    
    m = mp.Manager()
    out_q = m.Queue()
    nh = len(h.data)

    pool = mp.Pool(processes=ncore)
    for i in range(nh):
        pool.apply(mk_gal, args=(h.data[i], out_q, s.info, i,  work_dir))

    pool.close()
    pool.join()

    time.sleep(5)
    print("----------Done---------")
    print(" Took", time.time() - t1, "Seconds")

    dictout=[]
    try:
        f = open( work_dir + 'galaxies' + snout + '.txt', 'w')
    except:
        print("No filename is given.\n ")
    
    f.write(" #      ID        x          y       z[Mpc]       vx      vy     vz[km/s]")
    f.write("    Reff[kpc]      Mstar    Mgas[Msun]\n")    
    
    for i in range(nh):
        dd = out_q.get(timeout=2)
        f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(i,
                dd['id'],dd['xc'],dd['yc'],dd['zc']))
        f.write("  {:.3f}  {:.3f}  {:.3f}".format(dd['vx'],dd['vy'],dd['vz']))
        f.write("  {:.6f}  {:.0f}  {:.0f}     \n".format(dd['rgal'],dd['mstar'], dd['mgas']))
        dictout.append(dd)

    f.close()    
    
    catalog = pd.DataFrame(dictout).to_records()
    
    with open(fcat, 'wb') as f:
        pickle.dump(catalog, f)        
    
    print("\n", nout, "Done \n")
