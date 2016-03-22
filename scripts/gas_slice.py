# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 19:57:02 2015

@author: hoseung
"""

import load
import numpy as np
import utils.sampling as smp
import tree.halomodule as hmo 

#wdir = './'
#wdir = '/home/hoseung/Work/data/AGN2/'
wdir = '/home/hoseung/Work/data/Aquarius/'

rscale = 5.0
lmax = 10
#nout = int(input("NOUT? \n"))
nout=193

s = load.sim.Sim()
s.setup(nout, wdir)    

hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=s.info)
hh.load() 
#%%
i_center = np.where(hh.data['np'] == max(hh.data['np']))[0]

import draw
import utils.util
utils.util.reimport(draw)

#region = smp.set_region_multi(xc=hh.data.x[i_center],
#                              yc=hh.data.y[i_center],
#                              zc=hh.data.z[i_center],
#                              radius = hh.data.rvir * rscale)
region=smp.set_region()
zmi = region["ranges"][2][0]
zma = region["ranges"][2][1]

#n_slices = 100
n_slices = (zma - zmi) * s.info.pboxsize / 0.05 # one slice = 50kpc
dz = (zma - zmi) / n_slices
depth_in_cm = dz * s.info.pboxsize * 3.084e24 # 1Mpc = 3.084e24 cm
npix=2000
ptype = 'gas_den'


import matplotlib.pyplot as plt
import numpy as np
import utils.prettyplot as ptt
from matplotlib.colors import LogNorm

vmax = 1e-1
vmin = 1e-7

nticks = 5
dpi = 400
#cname ='gist_ncar'
cname ='gist_rainbow'
show=False
#%%
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    # By Jake VanderPlas
    # License: BSD-style

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

N_color = 14
import utils.util
utils.util.reimport(draw.pp)
column = False
#import pyximport; pyximport.install()
for i in range(int(np.ceil(n_slices))):
#for i in range(51, 52):
    region["ranges"][2]=[zmi + i * dz, zmi + (i+1) * dz]
    s.set_ranges(region["ranges"])
    #s.show_cpus()
    s.add_hydro()
    s.hydro.amr2cell(lmax=lmax)

    field = draw.pp.pp_cell(s.hydro.cell, npix, s.info, verbose=False, column=column)
    if column: 
        field /= depth_in_cm
    for j in range(2):
        save = wdir+str(i).zfill(3) + '_' + str(j) + 'hydro.png'
                        
        plt.ioff()
        fig = plt.figure()
        axes = fig.add_subplot(111)
        if j == 0:
            p = axes.imshow(field, cmap=discrete_cmap(N_color, cname),
                        norm=LogNorm(vmin=vmin, vmax=vmax))
        if j == 1:
            p = axes.imshow(field, cmap='gist_ncar',
                        norm=LogNorm(vmin=vmin, vmax=vmax))            
        xr = np.asarray(s.ranges[0])
        yr = np.asarray(s.ranges[1])
    
        x_old, x_new = ptt.tickInPboxInt(xr * s.info.pboxsize, field.shape[0], nticks=nticks)
        y_old, y_new = ptt.tickInPboxInt(yr * s.info.pboxsize, field.shape[1], nticks=nticks)
    
        fig.suptitle("Gas density averaged along the z-direction")
        plt.xticks(x_old * field.shape[0], x_new)
        plt.yticks(y_old * field.shape[1], y_new) # xticks, instead of set_xticks!
        #axes.set_xticks(x_old * field.shape[0], x_new)
    #    axes.set_yticks(y_old * field.shape[1], y_new)
    #            if zposition:
        annotation = 'z-range[Mpc]: {:.3f} - {:.3f}'.format(
                                region["ranges"][2][0]*s.info.pboxsize
                                ,region["ranges"][2][1]*s.info.pboxsize)
        axes.text(0.05, 1.02, annotation, transform = axes.transAxes, ha='left', fontsize=12)
    
        cbar = plt.colorbar(p)
        cbar.set_label(r'Hydrogen/$cm^{3}$', rotation=270, labelpad=10)
    
        if show:
            plt.show()
        if save:
            plt.savefig(save, dpi=dpi)
            plt.close()
    