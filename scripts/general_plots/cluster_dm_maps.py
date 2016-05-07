# -*- coding: utf-8 -*-
"""
DM density map of target clusters

Created on Tue Jun 23 03:46:41 2015

@author: hoseung
"""
import load
import tree
import utils.match as mtc
#%%
work_dir = '/home/hoseung/Work/data/DMO/'
nout=80

# load DMO 
s = load.sim.Sim(nout=nout, base=work_dir, dmo=True)
s.add_part(['dm id pos mass'])
s.set_ranges() # whole volume
s.part.load()
#%%
s.part.dm['m'] *= s.info.msun
#%%
# load halo
hall = tree.halomodule.Halo(nout=nout, base=work_dir, halofinder="HM", info=info)
hall.load()
ids = [32967,9782,5420,49919,36413,35663,32675,29176,24962,17891,10002,5427,
       4822,74010,49096,46814,41608,29195,29172,28928,28927,24954,14172,6098,
       49464,49367,41092,36370,36297,35861,35582,35401,28432,18212,18206,15233,
       5743,14168,28930,11402,26649,17545] # 43 target clusters

hind = mtc.match_list_ind(hall.data.id, ids)

h = tree.halomodule.Halo()
h.derive_from(hall, hind)

#%%
# Draw particle map to check.
#import zoom
import matplotlib.pyplot as plt
#fig, axs = plt.subplots(4,6)
#axs = axs.reshape(-1)

scale = 2.0
import draw
import utils.util 
utils.util.reimport(draw)
utils.util.reimport(draw.pp)
utils.util.reimport(draw.img_obj)

from matplotlib.colors import LogNorm

from mpl_toolkits.axes_grid1 import AxesGrid
plt.close()

"""
label_mode : 
1 = only one subplot
L = all left and all bottom

"""
from matplotlib.patches import Rectangle
import math

npix=400
vmin=1e13
vmax=1e20
#img_rad = 0.015 # 5Mpc/h on a side.
nrow=3
ncol=5
nclusters = len(h.data.id)
plotsperpage = nrow * ncol
npages = math.ceil(nclusters / plotsperpage)

# sort in decending mass order
ind=np.arange(nclusters)
ind = ind[np.argsort(h.data.mvir)[::-1]]

for page in range(0, npages):
    i_ini = page * (plotsperpage)
    i_fi = (page + 1) * (plotsperpage)

    fig = plt.figure(figsize=(10, 7))
    fig.subplots_adjust(left=0.01, right=0.99)
    grid = AxesGrid(fig, 111, # similar to subplot(142)
                        nrows_ncols = (nrow, ncol),
                        axes_pad = 0.05,
                        share_all=True,
                        label_mode = "1")
    """
                        ,cbar_location = "top",
                        cbar_mode="single",
                        cbar_size="3%"
                        ) 
    """
    for i in range(i_ini, min([i_fi, nclusters])):
        print('\n')
        print(page, i)
        img_rad = h.data.rvir[ind[i]] * 2
        h_region = smp.set_region(xc=h.data.x[ind[i]],
                                  yc=h.data.y[ind[i]],
                                  zc=h.data.z[ind[i]],
                                  radius = img_rad)
    #                            radius = img_rad) # +/- 3Mpc/h
        img = draw.pp.part2den(s.part.dm, info, region=h_region, npix=npix, verbose=True)
        
        im = grid[i-i_ini].imshow(img.data, norm=LogNorm(vmin=vmin, vmax=vmax),
                             interpolation="nearest", cmap='Greys')
    
        """
        # virial radius
        rad = npix * h.data.rvir[i] / img_rad
        circle = plt.Circle((npix/2,npix/2), radius=rad, fill=False)
        grid[i].add_patch(circle)   
        """
        # scalebar
        barlen = npix *0.005 / img_rad
        bar = Rectangle((0.1*npix, 0.9*npix), barlen, 12, facecolor="grey")
        grid[i-i_ini].add_patch(bar)

#        grid[i-i_ini].set_title(str(h.data.id[ind[i]]), fontsize=10)
        ann = grid[i-i_ini].annotate(str(h.data.id[ind[i]]),
                xy=(0.7*npix, 40),  xycoords='data',
                xytext=(0, 0), textcoords='offset points',
                size=10,
                bbox=dict(fc="white", ec="none")                
                )

        grid.axes_llc.set_xticks([])
        grid.axes_llc.set_yticks([])
        # Colorbar    
#        grid.cbar_axes[i-i_ini].colorbar(im)
#        grid.cbar_axes[i-i_ini].axis["top"].toggle(ticklabels=True,label=True)
#        grid.axes_llc.set_axis_off()
#        grid.axes_llc.set_xticks([npix*(1/2 - 2/5), npix/2 , npix*(1/2 + 2/5)])
#        grid.axes_llc.set_yticks([npix*(1/2 - 2/5), npix/2 , npix*(1/2 + 2/5)])
#        grid.axes_llc.set_xticklabels([str(i) for i in range(-2, 3, 2)])
#        grid.axes_llc.set_yticklabels([str(i) for i in range(-2, 3, 2)])
#        grid.axes_llc.set_xlabel('Mpc/h')
#        grid.axes_llc.set_ylabel('Mpc/h')
        
    fig.savefig(work_dir + 'clusters_DM_' + str(page) +'.png', dpi=160)
    plt.close()

plt.show()
#plt.close()


#%%
def grid_with_single_cbar(fig):
    """
    A grid of 2x2 images with a single colorbar
    """
    grid = AxesGrid(fig, 142, # similar to subplot(142)
                    nrows_ncols = (2, 2),
                    axes_pad = 0.0,
                    share_all=True,
                    label_mode = "L",
                    cbar_location = "top",
                    cbar_mode="single",
                    )

    Z, extent = get_demo_image()
    for i in range(4):
        im = grid[i].imshow(Z, extent=extent, interpolation="nearest")
    #plt.colorbar(im, cax = grid.cbar_axes[0])
    grid.cbar_axes[0].colorbar(im)

    for cax in grid.cbar_axes:
        cax.toggle_label(False)

    # This affects all axes as share_all = True.
    grid.axes_llc.set_xticks([-2, 0, 2])
    grid.axes_llc.set_yticks([-2, 0, 2])
