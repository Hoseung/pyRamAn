# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 12:48:45 2015

@author: hoseung
""" 


def mk_gal(halodata, out_q, info, i, final_gal,
            save=False, rscale=0.3, verbose=False, galaxy_plot_dir='./',
            rscale_lambda=2.0, npix_lambda=50, npix=400, galaxy_plot=False,
            method_com=2):
    """
    Direct plot,
    Create galaxy, 
    Calculate lambda_r (using Cappellari 2003)
    Draw ed map of galaxy.
    
    """
    from galaxy import galaxy
    import utils.sampling as smp
    import draw
    import matplotlib.pyplot as plt    
#    print("IDs:", id(star), id(dm), id(cell))

    gal_out = {"id":0, "xc":0.0, "yc":0.0, "zc":0.0,
               "vx":0.0, "vy":0.0, "vz":0.0,
               "mstar":0.0, "nstar":0.0, "mgas":0.0,
               "lambda_arr":[], "lambda_r":0, "rgal":0, "final":final_gal,
                "rhalo":halodata['rvir'], "boxtokpc":info.pboxsize*1000}
    # Direct plot ---------------------------------------------------------
    region = smp.set_region(xc=halodata['x'],
                            yc=halodata['y'],
                            zc=halodata['z'],
                            radius = halodata['rvir'])
    if galaxy_plot:
        extent = (0, npix, 0, npix)        
        star_map = draw.pp.den2d(star['x'],star['y'],star['z'],star['m'], npix,
                                 region=region, cic=True, norm_integer=False)
        if star_map is not False:
            ls = np.zeros((npix,npix))
            ii = star_map > 0
            ls[ii] = np.log10(star_map[ii]) # Stellar map HAS empty pixels.
            ls[star_map <= 0] = np.floor(ls[ii].min())
            plt.imshow(ls, cmap="CMRmap", interpolation='nearest', extent=extent)
        
        # One of two should be transposed.
        # But which one?
#        gas_map = draw.pp.pp_cell(cell, npix, info, region=region, verbose=False)
#        im2 = plt.imshow(np.transpose(np.log10(gas_map)), cmap="CMRmap", alpha=.5, interpolation='bilinear', extent=extent)
    
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
        
        plt.savefig(galaxy_plot_dir+"2dmap_"+str(halodata['id']).zfill(5)+'.png', dpi=144)
        plt.close()

    #Create galaxy ---------------------------------------------------------
    gal = galaxy.Galaxy(halodata, radius_method='eff', info=info)
    is_gal = gal.mk_gal(star=star, dm=dm, cell=cell,
               rscale=rscale, verbose=verbose, method_com=method_com)
    #-----------------------------------------------------------------------    
    print(gal.id, "IS_GAL",is_gal)
    if not is_gal:
        print(gal.id, " Not a good galaxy")
        out_q.put(gal_out)
    else:
        # Save to catalog -------------------------------------------------------
#        print("Good galaxy, R eff:", gal.reff)
        gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    
        # Calculate lambda_r ---------------------------------------------------

        gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                             + "_" + str(final_gal).zfill(5) + "_"  \
                             + str(gal.id) + ".png", ioff=True)
#       gal.save_gal(base=wdir)

        # save in a dict.
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
    return gal


def plot_lambda(catalog, i_early, i_late, i_bad, fn_out='./'):
    import matplotlib.pyplot as plt
    plt.ioff()
    f = plt.figure()
    ax = f.add_subplot(111)
    #for i, val in enumerate(lambdar_arr):
    for i in i_early:
        a = np.asarray(catalog['lambda_arr'][i])
        ax.plot(a, 'r-', alpha=0.5) # Red = Early
    for i in i_late:
        ax.plot(catalog['lambda_arr'][i], 'b-', alpha=0.3) # Red = Early
    
    #plt.xlabel() # in the unit of Reff
    ax.set_title(r"$\lambda _{R}$") 
    ax.set_ylabel(r"$\lambda _{R}$") 
    ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
    ax.set_xlim(right=9)
    ax.set_xticks([0, 4.5, 9])
    ax.set_xticklabels(["0", "0.5", "1"])
    plt.savefig(fn_out)
    plt.close()    

"""
The processing pool needs to be instantiated in the main 
thread of execution. 
"""
import load
from tree import tmtree
import numpy as np
import utils.sampling as smp
import tree.halomodule as hmo 
import os

#wdir = input("Working directory \n")
wdir = '/home/hoseung/Work/data/05427/'

#wdir = './'
nout_fi = 187 # fixed. Tree depends on this value. 
final_gal = int(input("idgal: \n"))
nout_ini = input("First snapshot: (default = 37 , z=3) \n")
nout_end = input("Last snapshot: (default = 187, z=0) \n")
#final_gal = 1392
#nout_ini=130
#nout_fi=187
if final_gal == "":
    import sys
    sys.exit(0)

if nout_ini == "":
    nout_ini = 37
else:
    nout_ini = int(nout_ini)

if nout_end == "":
    nout_end = 187
else:
    nout_end = int(nout_end)

#nout_end = 131


#----------------------------------------------------------------------
# 27 : z=4;  37 : z=3;  20 : ~ z=5
nouts = range(nout_fi, nout_ini -1, -1) 
#----------------------------------------------------------------------

mstar_min_plot = 1e9
mk_gal_rscale = 1.0 # unit of Rvir,halo
rscale = 1.5
# should show fly-bys..

npix=800
rscale_lambda = 3.0 # Reff unit.
npix_lambda = int(10 * rscale_lambda)
lmax = 19
ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]

# optional parameters ----------------------------------------------------
lambda_plot = False 

## halo -----------------------------------------------------------
m_halo_min = 1e10 # minimum halo mass above which galaxies are searched for. 
dir_out = wdir + 'catalog/'

if is_gal:
    # Galaxy tree
    tree_path = 'GalaxyMaker/Trees/'
else:
    # halo tree
    tree_path = 'halo/Trees/'
    
try:    
    alltrees = pickle.load(open(wdir + tree_path + "extended_tree.pickle", "rb" ))
    print("Loaded an extended tree")
except:
    alltrees = treemodule.CTree()
    alltrees.load(filename= wdir + tree_path + 'tree_0_0_0.dat')
    # Fix nout -----------------------------------------------------
    nout_max = alltrees.data['nout'].max()
    alltrees.data['nout'] += nout_fi - nout_max
    print("------ NOUT fixed")
    alltrees.data = ctu.augment_tree(alltrees.data, wdir, is_gal=is_gal)
    print("------ tree data extended")
    pickle.dump(alltrees, open(wdir + tree_path + "extended_tree.pickle", "wb" ))

tt = alltrees.data
tt_final = tt[tt['nout'] == nout_fi]

info = load.info.Info(nout=nout_fi, base=wdir, load=True)
# Get the tree of the target halo(s)
# No need to check if the halo is near the clsuter or not. 
prg_idx, prg_id = tmtree.get_main_prg(tt, final_gal)
#%%
#import galaxy
import utils.util
utils.util.reimport(load)
utils.util.reimport(load.utils)
try:
    if not os.path.isdir(dir_out):
        os.mkdir(dir_out)
    f = open(dir_out + 'galaxies' +str(final_gal) + '.txt', 'w')
except:
    print("No filename is given.\n ")

f.write(" Final ID:" + str(final_gal) + "\n")
f.write(" nout     ID        x          y       z[Mpc]       vx      vy     vz[km/s]")
f.write("    Reff[kpc]     Mstar    Mgas[Msun]  Rhalo[kpc]  boxtokpc  \n")    

print("!!!!!!!!!!!!!!!!!!")
import queue
out_q = queue.Queue()
dictout=[]
dd = None
galaxy_plot_dir = wdir + 'galaxy_' + str(final_gal).zfill(5) + '/'
if os.path.isdir(galaxy_plot_dir) is False:
    os.mkdir(galaxy_plot_dir)


#%%
r_mod = 1.0
for inout, nout in enumerate(nouts):
    if nout > nout_end:
        continue
    print(inout, nout, nout_end)
    snout = str(nout)
    
    info = load.info.Info(nout=nout, base=wdir, load=True)
           
    # Load all halo
    hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info, load=True, is_gal=True)
#    hind = match.match_list_ind(hh.data['id'], halo_ok[igal,inout])   
    h = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info)
    h.derive_from(hh, hh.data['id'] == prg_id[0,inout])

    if dd is not None and 'rhalo' in dd.keys():
        rhal_min_possible = max([h.data['rvir'], 0.9 * r_mod * dd['rhalo']])
        r_mod = rhal_min_possible/(h.data['rvir'] )
        print(" Rscale multiplied by ", r_mod)
    
    r_mod = 1.0
#    print(h.data)
    region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z,
                                  radius = h.data.rvir * rscale * r_mod)

    s = load.sim.Sim(nout=nout, base=wdir, ranges=region["ranges"], setup=True)
#    s.show_cpus()
    s.add_part(ptypes, load=True, fortran=True)
    assert s.part.nstar > 0, "Not enough stellar particles in given cpus"
        
    s.add_hydro(load=True, lmax=lmax)

    star = s.part.star
    dm = s.part.dm
    cell = s.hydro.cell

    keywords = dict(galaxy_plot_dir=galaxy_plot_dir,
                    rscale = mk_gal_rscale,
                    verbose=True, rscale_lambda=rscale_lambda,
                    npix_lambda=npix_lambda, galaxy_plot = True,
                    method_com=2)
    # To keep the consistency, use queue to get output data. 

    gal = mk_gal(h.data[0], out_q, info, 0, final_gal, **keywords)

    tmp =  out_q.get()
    if tmp['id'] == 0:
        continue
    dd = tmp
    f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(nout,
            dd['id'],dd['xc'],dd['yc'],dd['zc']))
    f.write("  {:.3f}  {:.3f}  {:.3f}".format(dd['vx'],dd['vy'],dd['vz']))
    f.write("  {:.6f}  {:.0f}  {:.0f}".format(dd['rgal'],dd['mstar'], dd['mgas']))
    f.write("  {:.5f}  {:.5f}    \n".format(dd['rhalo'],dd['boxtokpc']))
    dictout.append(dd)
    print("  ")

print("----------Done---------")
f.close()    
print("Text file written")

import pandas as pd
catalog = pd.DataFrame(dictout).to_records()
# pickling also works!
# The problem was memory being locked!
import pickle
fcat = dir_out +"catalog" + str(nout_fi) + "_" + str(final_gal) + ".pickle"
with open(fcat, 'wb') as f:
    pickle.dump(catalog, f)  


#%%

if lambda_plot:
    i_early = np.where(catalog['mstar'] > mstar_min_plot)[0]
    i_late = []
    i_bad = np.where(catalog.id == 0)[0]    
    plot_lambda(catalog, i_early, i_late, i_bad,
                fn_out = snout + "_lambdar_disk.png")
