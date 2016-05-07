# -*- coding: utf-8 -*-
"""
Plots merger history

First used for SAMI project.
This is also useful in showing merging cluster's merger history.

Created on Thu Jun  4 22:58:03 2015
@author: hoseung
"""
# Load a tree 
from astropy.io import fits
from astropy.table import Table

import load
import tree
import numpy as np
import utils.sampling as smp
from utils import util
from tree import tmtree
import tree.draw_merger_tree as dmt
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

wdir = '/home/hoseung/Work/data/05427/'
#data = fits.getdata(wdir + "halo/TMtree.fits", 1)
data = fits.getdata(wdir + "GalaxyMaker/TMtree.fits", 1)
tt = Table(data)

#%% 

nout_ini = 37
nout_fi = 187

info = load.info.Info(nout=nout_fi, base=wdir)

frefine= 'refine_params.txt'
fnml = 'cosmo_200.nml'

ptypes=["star id pos mass vel", "dm id pos mass vel"]

# Load all halo
#hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="RS", info=info)
hall = tree.halomodule.Halo(nout=nout_fi, base=wdir, halofinder="HM", info=info,
                            is_gal=True, load=True)

# convert to code unit. - done by default
#hall.normalize()

# subset of halos ONLY inside zoom-in region
i_center = np.where(hall.data['np'] == max(hall.data['np']))[0]

h_ind = smp.extract_halos_within(hall.data, i_center, scale=2.0)

#%%
h = tree.halomodule.Halo()
h.derive_from(hall, h_ind)
region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir)

# filter by volume

hind = np.where( (hall.data.x > region["xr"][0]) & (hall.data.x < region["xr"][1]) &
                (hall.data.y > region["yr"][0]) & (hall.data.y < region["yr"][1]) &
                (hall.data.z > region["zr"][0]) & (hall.data.z < region["zr"][1]) &
                (hall.data.mvir > 1e11))[0]

                
h.derive_from(hall, hind)
#%%
# pick galaxies to track. (one at first)

# check_tree_complete depends on tmtree.get_main_prg.
# I see no reason this should be halofinder-specific. 
# make it general! 
#inout = np.where(tt[:]["NOUT"] == 0)[0]
#halo_list = h.data["id"][1:]
#h_ind_ok, halo_ok = tmtree.check_tree_complete(tt, 0, nout_fi - nout_ini, halo_list)

#thehalo=halo_ok[0][50] # 1993 , IDX = 1054465
#%%
# fig 1. number of mergers 

#tfinal = tt[inout]
#ihal = np.where(tfinal['HALNUM'] == thehalo)[0]
#print(tfinal['IDX'][ihal])



#%%
# fig 2. merger tree 
#ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))

#for galid in h.data.id[0:10]:
ttt = tt[tt["NOUT"] == 0] # final galaxies
for galid in ttt["HALNUM"]:
    print(galid)
    fig = plt.figure(figsize=[6,6])
    plt.ioff()
    ax = fig.add_subplot(111)
    #galid = 6033
    sidgal = str(galid).zfill(5)
    
    nouts = np.unique(range(nout_fi - nout_ini + 2))
    zreds = np.unique(tt["Z"])[:len(nouts)]
    zreds = ["%.2f" % i for i in zreds]
    #print(zreds)
    idx0 = dmt.get_idx(tt, galid, 0)
    dmt.recursive_tree(idx0, tt, 123, ax, 0, 0, 1, mass_unit=1e9)
    
    # y axis label (redshift)
    ax.set_ylabel("Redshift")
    ax.set_xlim([-1.1,2])
    ax.set_ylim([0,151])
    plt.yticks(nouts[1:151:10], zreds[1:151:10])
    ax.set_title(sidgal)
    #fig.show()
    plt.savefig(wdir + "mergertrees/" + sidgal + '.png')
    plt.close()

