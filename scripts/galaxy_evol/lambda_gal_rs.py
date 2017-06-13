# coding: utf-8

def mk_gal(halodata, out_q, info, final_gal,
            save=False, rscale=0.3, verbose=False, galaxy_plot_dir='./',
            rscale_lambda=2.0, npix_lambda=50, npix=400, galaxy_plot=False):
    """
    Direct plot,
    Create galaxy, 
    Calculate lambda_r (using Cappellari 2003)
    Draw ed map of galaxy.
    
    """
    from galaxymodule import galaxy
    import utils.sampling as smp
    import draw
    import matplotlib.pyplot as plt    
#    print("IDs:", id(star), id(dm), id(cell))

    gal_out = {"id":0, "xc":0.0, "yc":0.0, "zc":0.0,
               "vx":0.0, "vy":0.0, "vz":0.0,
               "mstar":0.0, "nstar":0.0, "mgas":0.0,
               "lambda_arr":[], "lambda_r":0, "rgal":0, "final":final_gal}
    # Direct plot ---------------------------------------------------------
    region = smp.set_region(xc=halodata['x'][0],
                            yc=halodata['y'][0],
                            zc=halodata['z'][0],
                            radius = halodata['rvir'][0])
    # Without trailing [0], region becomes list of list of arrays.
        # This causes error in x,y range setting in the pp_cell function.                            
    if galaxy_plot:
        print("Plotting region")
        extent = (0, npix, 0, npix)        
        star_map = draw.pp.den2d(star['x'],star['y'],star['z'],star['m'], npix,
                                 region=region, cic=True, norm_integer=False)
        if star_map is not False:
            ls = np.zeros((npix,npix))
            ii = star_map > 0
            ls[ii] = np.log10(star_map[ii]) # Stellar map HAS empty pixels.
            ls[star_map <= 0] = np.floor(ls.min())
            plt.ioff()
            im1 = plt.imshow(ls, cmap="gray",
                             interpolation='nearest',
                             extent=extent)       
        
        # One of two should be transposed.
        # But which one?

        gas_map = draw.pp.pp_cell(cell, npix, info,
                                  region=region,
                                  verbose=False)
        im2 = plt.imshow(np.transpose(np.log10(gas_map)),
                         cmap="CMRmap",
                         alpha=.5,
                         interpolation='bilinear',
                         extent=extent)
    
        rgal = region['radius'] * s.info.pboxsize * 1000

        ax = plt.gca()
        ax.set_xlabel("position [kpc]")
        ax.set_xticks(np.linspace(0,npix,5))
        xticks = ["{:.2f}".format(x) \
                    for x in np.linspace(-rgal, rgal, num=5)]
        ax.set_xticklabels(xticks)
        ax.set_ylabel("position [kpc]")
        ax.set_yticks(np.linspace(0,npix,5))
        yticks = ["{:.2f}".format(y) \
                    for y in np.linspace(-rgal, rgal, num=5)]
        ax.set_yticklabels(yticks)
        
        plt.savefig(galaxy_plot_dir + "2dmap_" + 
                    str(halodata['id'][0]).zfill(5)+'.png', dpi=144)
        plt.close()

    #Create galaxy ---------------------------------------------------------
    gal = galaxy.Galaxy(halodata, radius_method='simple', info=info)
    is_gal = gal.mk_gal(star=star, dm=dm, cell=cell,
               rscale=rscale, verbose=verbose)
    #-----------------------------------------------------------------------    
    print(gal.id, "IS_GAL",is_gal)
    if not is_gal:
        print(gal.id, " Not a good galaxy")
        out_q.put(gal_out)
    else:
        # Save to catalog -------------------------------------------------------
        gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    
        # Calculate lambda_r ---------------------------------------------------

        gal.plot_gal(save_dir = galaxy_plot_dir, ioff=True)
    #            gal.save_gal(base=wdir)

        # Instead of galaxy class, save them in a dict.    
        gal_out['mstar'] = gal.mstar
        gal_out['mgas'] = gal.mgas
        gal_out['nstar'] = gal.nstar
        gal_out['id'] = gal.id
        gal_out['xc'] = gal.xc * info.pboxsize
        gal_out['yc'] = gal.yc * info.pboxsize
        gal_out['zc'] = gal.zc * info.pboxsize
        gal_out['vx'] = gal.vxc * info.kms
        gal_out['vy'] = gal.vyc * info.kms
        gal_out['vz'] = gal.vzc * info.kms        
        gal_out['lambda_arr'] = gal.lambda_arr
        gal_out['lambda_r'] = gal.lambda_r
        gal_out['rgal'] = gal.reff# * info.pboxsize * 1000.0 # in kpc
        out_q.put(gal_out)
    print("mk_gal done")
    return gal_out
    
#%%
"""
The processing pool needs to be instantiated in the main 
thread of execution. 
"""
import multiprocessing as mp    
import load
from tree import tmtree
import numpy as np
import utils.sampling as smp
import tree.halomodule as hmo 
from utils import match
import os

#wdir = input("Working directory \n")
#nout = int(input("nout? \n"))
wdir = './'
wdir = '/home/hoseung/Work/data/AGN2/'
ncore = 1
dir_out = wdir + 'catalog/'
ftree = wdir + 'rhalo/tree.pickle'


#----------------------------------------------------------------------
# 27 : z=4;  37 : z=3;  20 : ~ z=5
nout_ini = 37
nout_fi = 132
nouts = range(nout_fi, nout_ini -1, -1) 
nnouts = nout_fi - nout_ini + 1
#----------------------------------------------------------------------

rscale = 0.8
r_cluster_scale = 2.0 # maximum radius inside which galaxies are searched for
npix=800
rscale_lambda = 3.0
npix_lambda = int(15 * rscale_lambda)
lmax = 19
ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]

## halo part -----------------------------------------------------------



# optional parameters ----------------------------------------------------

# nout_halo = 122 == nout 10, nout_halo = 0   == nout 132
TreeMaker = False
ConsistentTree = True
if TreeMaker:
    tt = tmtree.load(work_dir=wdir, filename="halo/TMtree.fits")
    tfin = tt[np.where(tt['NOUT'] == 0)]
    tini = tt[np.where(tt['NOUT'] == nout_fi - nout_ini)]

if ConsistentTree:
    import tree.treemodule as tmo
    from utils.util import reimport
    ftree = wdir + 'rhalo/tree.pickle'
    reimport(tmo)
    tt = tmo.CTree(filename=ftree)

#%%
import tree.treeutils as tru
reimport(tru)

halofinder = "RS"    

roots = np.unique(tt.data['tree_root_id'])
igal = 0 # for only one galaxy.
#target_halo = roots[227153]
target_halo = [1630435]

#h_ind_ok, halo_ok = tmtree.check_tree_complete(tt, 0, nout_fi - nout_ini, target_halo)
h_ind_ok, halo_ok = tru.check_tree_complete(tt, target_halo,
                                            nout_ini = 0, nout_fi =81)

# Be ware taht halo_ok list is unique tree ID, 
# whereas hh.data['id'] == tt.data['Orig_halo_id'].
#%%
halo_ok_org=np.zeros(len(halo_ok[0]), dtype=int)
for i, this_hal in enumerate(halo_ok[0]):
    i_now = np.where(tt.data['id'] == this_hal)[0]
    halo_ok_org[i] = tt.data['Orig_halo_id'][i_now]

#%%
#h_ind_ok, halo_ok = tru.check_tree_complete(tt.data, 0, nout_fi - nout_ini, target_halo)
print(len(halo_ok), "halo(s) left")
 
gal_history = []
sgalid = str(target_halo[0]).zfill(5)
fcat = dir_out + sgalid + "catalog" + ".pickle"
try:
    if not os.path.isdir(dir_out):
        os.mkdir(dir_out)
    f = open(dir_out + 'galaxies' + sgalid + '.txt', 'w')
    f.write(" nout      zred      ID        x          y")
    f.write("       z[Mpc]       vx      vy     vz[km/s]")
    f.wrtie("    Reff""[kpc]      Mstar    Mgas[Msun] ")
    f.write("\n")
except:
    print("No filename is given.\n ")

#for inout, nout in enumerate(nouts):
final_gal = target_halo[0]
for inout, nout in enumerate([132]):
    print(inout, nout)
    snout = str(nout)
   
    info = load.info.Info(nout=nout, base=wdir, load=True)
           
    # Load all halo
    if halofinder == "HM":
        hh = hmo.Halo(base=wdir, nout=nout,
                      halofinder=halofinder, info=info,
                      load=True)   
        # Find a halo with halo_ok.    
        hind = np.where(hh.data['id'] == halo_ok_org[inout])
        # assume halo_ok.shape = (1,nnout)
        h = hmo.Halo(base=wdir, nout=nout, halofinder=halofinder, info=info)
        h.derive_from(hh, hind)
        region = smp.set_region_multi(xc = h.data.x,
                                      yc = h.data.y,
                                      zc = h.data.z,
                                      radius = h.data.rvir * rscale)
    
    if  halofinder == "RS":
        h = hmo.Halo(base=wdir, nout=nout,
                      halofinder=halofinder, info=info,
                      load=True)
        #ihal = np.where(tt.data['id'] == halo_ok[igal, inout])
        ihal = np.where(h.data['id'] == halo_ok_org[inout])
        # tt.data['x,y,z'] in Mpc/h unit. 
        # tt.data['rvir'] in kpc/h unit. 
        # convert it to code unit.
#        region = smp.set_region_multi(\
#            xc = tt.data['x'][ihal]/199.632011,
#            yc = tt.data['y'][ihal]/199.632011,
#            zc = tt.data['z'][ihal]/199.632011,
#            radius = tt.data['rvir'][ihal]/199.632011/1000.0 * rscale)
        region = smp.set_region(xc = h.data['x'][ihal][0],
                                yc = h.data['y'][ihal][0],
                                zc = h.data['z'][ihal][0],
                                radius = h.data['rvir'][ihal][0] * rscale)
        # Without trailing [0], region becomes list of list of arrays.
        # This causes error in x,y range setting in the pp_cell function.

    s = load.sim.Sim(nout=nout, base=wdir, ranges=region["ranges"], setup=True)
#    s.show_cpus()
    s.add_part(ptypes, load=True, fortran=True)
    assert s.part.nstar > 0, "Not enough stellar particles in given cpus"
        
    s.add_hydro(load=True, lmax=lmax)

    star = s.part.star
    dm = s.part.dm
    cell = s.hydro.cell

#   Multiprocessing -----------------------------------------------------------
    m = mp.Manager()
    out_q = m.Queue()
#    nh = len(h.data)

    keywords = dict(galaxy_plot_dir=dir_out,
                    verbose=False, rscale_lambda=rscale_lambda,
                    npix_lambda=npix_lambda, galaxy_plot = True)
    print("  ")
    print("mk_gal start")
    mk_gal(h.data[ihal], out_q, s.info, final_gal, **keywords)
#    print("Looking for galaxies inside {} halos".format(nh))
#    pool = mp.Pool(processes=ncore)
#    for i in range(nh):
#        pool.apply_async(mk_gal, args=(h.data[i], out_q,
#                                       s.info, i, halo_ok[i,0]), kwds=keywords)
    # Pass optional parameters in dict!s
#    pool.close()
#    pool.join()    
    print("----------Done---------")

    try:
        tmp =  out_q.get(timeout=1)
        if tmp['id'] == 0:
            continue
        dd = tmp
        f.write("{:<3}   {:.5f}    {:<4}   {:.5f}  {:.5f}  {:.5f}".format(
                nout, info.zred, dd['id'],dd['xc'],dd['yc'],dd['zc']))
        f.write("  {:.3f}  {:.3f}  {:.3f}".format(dd['vx'],dd['vy'],dd['vz']))
        f.write("  {:.6f}  {:.0f}  {:.0f}   \n".format(dd['rgal'],dd['mstar'],
                dd['mgas']))
        gal_history.append(dd)
    except:
        continue


#%%-----------------------
f.close()
print("Text file written")

import pandas as pd
catalog = pd.DataFrame(gal_history).to_records()    
# pickling also works!
import pickle
with open(fcat, 'wb') as f:
    pickle.dump(catalog, f)  
    
