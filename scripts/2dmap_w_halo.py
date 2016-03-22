# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:11:20 2015

@author: hoseung
"""



#sim_dir = ['C29195', 'C01605_hi/', 'C04466/', 'G2/kisti'][3]
#base_dir = path.join(base, '', sim_dir, '')
base_dir = '/home/hoseung/Work/data/05427/'

nout=186

import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(111)


import pickle
snout = str(nout).zfill(3)
fin = base_dir + snout + 'map.pickle'
with open(fin, 'rb') as f:
    img = pickle.load(f)

ptimg = img.ptden2d
fout = base_dir + snout + "dmmap_" + ptimg.proj + ".png"
img.ptden2d.plot_2d_den(save=False, show=False, vmin=1e13, vmax=1e20, dpi=200, axes=ax1)
#%%

import tree
import load
import numpy as  np
import utils.sampling as smp

#s = load.sim.Sim(nout, base_dir)
info = load.info.Info(nout=nout, base=base_dir, load=True)
hall = tree.halomodule.Halo(nout=nout, base=base_dir, halofinder="HM", info=info)
hall.load()

i_center = np.where(hall.data['np'] == max(hall.data['np']))
h = tree.halomodule.Halo()
h.derive_from(hall, [i_center]) 

region = smp.set_region(xc=h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * 2)      

#%%
from draw import pp
npix = 400

ind = np.where(hall.data.mvir > 5e11)
h_sub = tree.halomodule.Halo()
h_sub.derive_from(hall, ind) 
#x = hall.data.x#[ind]
#y = hall.data.y#[ind]
#r = hall.data.rvir#[ind]
#pp.circle_scatter(ax1, x*npix, y*npix, r*30, facecolor='none', edgecolor='b', label='555')

#ax1.set_xlim(right=npix).
#ax1.set_ylim(top=npix)
pp.pp_halo(h_sub, npix, region=img.ptden2d.region, axes=ax1, rscale=40, name=True)

plt.show()



#%%
