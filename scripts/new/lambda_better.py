# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 05:06:15 2015

@author: hoseung
"""

import time
import numpy as np

def extract_data(halo, rscale=1.25, radius='rvir'):
    xc_tmp0 = halo['x']
    yc_tmp0 = halo['y']
    zc_tmp0 = halo['z']
       
    #rr_tmp0 = min([halo[radius] * rscale, 0.0002]) 
    rr_tmp0 = max([halo[radius] * rscale, 0.0001])
    # arbitrary! < 20kpc/h
    
    # When merger occurs, larger radius is likely to include 
    # companion galaxy resulting center to be in the middle of nowhere.
    # If you want a larger galaxy, # increase rgal_tmp instead. 
    #        
    # xx is easier to search for than x.

    if star_all is not None:
        ind_s = np.where((star_all['x'] - xc_tmp0)**2 + (star_all['y'] - yc_tmp0)**2 
                        + (star_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    if dm_all is not None:
        ind_d = np.where((dm_all['x'] - xc_tmp0)**2 + (dm_all['y'] - yc_tmp0)**2 
                        + (dm_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    if cell_all is not None:
        ind_c = np.where((cell_all['x'] - xc_tmp0)**2 + (cell_all['y'] - yc_tmp0)**2 
                        + (cell_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    else:
        return star_all[ind_s], dm_all[ind_d], None
        
    return star_all[ind_s], dm_all[ind_d], cell_all[ind_c]



def mk_gal(halodata, info, i, idx,
           save=False, rscale=1.5, verbose=False, galaxy_plot_dir='./',
           rscale_lambda=3.0, npix_lambda=50, npix=400, galaxy_plot=False,
           method_com=2, mstar_min=5e9):
    
    from galaxy.galaxy2 import Galaxy

    rscale_extract = 1.1
    star, dm, cell = extract_data(h.data[i], rscale=rscale_extract, radius='rvir')
    if mstar_min > 0.:
        if sum(star['m']) * info.msun < mstar_min:
            print("(1)Not enough stars: {:.2f} Msun".format(sum(star['m']) * info.msun))
            print("Aborting... \n")
            print(" Not a good galaxy")
            out_q.put(gal_out)
            return
               
    #Create galaxy ----------------------------------------------
    gal = Galaxy(halodata, radius_method='eff', info=info)
    #print(i, time.time() - t, "seconds ---2")
    good_gal = gal.mk_gal_from_gal(star, dm, cell,
                        mstar_min=mstar_min,
               rscale=min([rscale,rscale_extract]), verbose=verbose, method_com=method_com)

    
    gal.cal_lambda_better(npix_lambda=5, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    

    gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                         + "_" + str(idx).zfill(7) + "_"  \
                         + str(gal.id) + ".png", ioff=True)
                             
    return gal
    

"""
The processing pool needs to be instantiated in the main 
thread of execution. 
"""
import load
import matplotlib.pyplot as plt
from tree import tmtree, treemodule
import tree.ctutils as ctu
import utils.sampling as smp
import tree.halomodule as hmo 
import os
import pickle
hydro = True
is_gal = True

wdir = '/home/hoseung/Work/data/05427/'
nout = 187

#----------------------------------------------------------------------
mstar_min = 5e9
# Only galaxies above this stellar mass at the final snapshot are considered.

mk_gal_rscale = 2.0 # unit of Rvir,galaxy
rscale = 1.5
npix=800
rscale_lambda = 3.0 # Reff unit.
npix_lambda = int(10 * rscale_lambda)
lmax = 19
ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]

## halo part -----------------------------------------------------------
m_halo_min = 5e9 # minimum galaxy mass above which galaxies are searched for. 
dir_out = wdir + 'catalog_GM/'

# optional parameters ----------------------------------------------------
do_galaxy_plot=False
lambda_plot = False

info = load.info.Info(nout=nout, base=wdir, load=True)

halo_ok = 76

# Do I really need halo bricks?
hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info, load=True, is_gal=is_gal)
hind = hh.data['id'] == halo_ok
h = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info, is_gal=is_gal)
h.derive_from(hh, hind)
halodata = h.data
region = smp.set_region(xc=halodata['x'],
                    yc=halodata['y'],
                    zc=halodata['z'],
                    radius = halodata['rvir'])


s = load.sim.Sim(nout=nout, base=wdir, setup=True, region = region)
#s.set_ranges(region = region)
s.add_part(ptypes, load=True, fortran=True)

if hydro:
    s.add_hydro(load=True, lmax=lmax)
    cell_all = s.hydro.cell
else:
    cell_all = None

star_all = s.part.star
dm_all = s.part.dm

keywords = dict(galaxy_plot_dir='./',
            rscale = mk_gal_rscale,
            verbose=False, rscale_lambda=rscale_lambda,
            npix_lambda=npix_lambda, galaxy_plot = do_galaxy_plot,
            mstar_min=mstar_min)


from queue import Queue
out_q = Queue()

gal = mk_gal(h.data, s.info, 0, 12345, **keywords)

print("----------Done---------")

from draw import pp
from matplotlib.colors import LogNorm
# npix = round(gal.nstar**(1/3) * 3) # to maintain roughly same pixel density. 
# Or, constant physical scale

#rscale = 3

#npix = round(gal.reff * rscale) * 4 # 5 * reff = Rgal in most case, 4 pixels in 1 kpc.
#region = smp.set_region(xc=0, yc=0, zc=0, radius = gal.reff * rscale)  
#data = pp.den2d(gal.star['x'], gal.star['y'], gal.star['z'], gal.star['m'], \
#              npix, region=region, cic=True, norm_integer=True)
              
              
#gal.cal_lambda_better()
plt.show()
#%%
gal.plot_gal()