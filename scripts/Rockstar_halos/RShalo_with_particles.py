
# coding: utf-8

# In[2]:

def plot_halo_region(region, info, star, cell=None,
                     npix=400,
                     fn_save='region_map.png',
                     show=True):
    import draw  
    extent = (0, npix, 0, npix)

    star_map = draw.pp.den2d(star['x'],star['y'],star['z'],star['m'], npix,
                             region=region, cic=True, norm_integer=False)
    if star_map is not False:
        ls = np.zeros((npix,npix))
        ii = star_map > 0
        ls[ii] = np.log10(star_map[ii]) # Stellar map HAS empty pixels.
        ls[star_map <= 0] = np.floor(ls.min())
        im1 = plt.imshow(ls, cmap="gray",
                         interpolation='nearest',
                         extent=extent,
                         origin='lower')
    
    # One of two should be transposed.
    # But which one?
    if cell is not None:
        gas_map = draw.pp.pp_cell(cell, npix, info, region=region, verbose=False)
        im2 = plt.imshow(np.transpose(np.log10(gas_map)),
                         cmap="CMRmap",
                         alpha=.5,
                         interpolation='bilinear',
                         extent=extent,
                         origin='lower')

    rgal = region['radius'] * info.pboxsize * 1000

    ax = plt.gca()
    ax.set_xlabel("position [kpc]")
    ax.set_xticks(np.linspace(0,npix,5))
    xticks = ["{:.2f}".format(x)                 for x in np.linspace(-rgal, rgal, num=5)]
    ax.set_xticklabels(xticks)
    ax.set_ylabel("position [kpc]")
    ax.set_yticks(np.linspace(0,npix,5))
    yticks = ["{:.2f}".format(y)                 for y in np.linspace(-rgal, rgal, num=5)]
    ax.set_yticklabels(yticks)
    
    plt.savefig(fn_save, dpi=144)
    if show:
        plt.show()
    else:
        return


# In[3]:

import numpy as np
import load
import tree
import matplotlib.pyplot as plt
import utils.sampling as smp
from draw import pp


# In[4]:

wdir = '/media/hoseung/btrfs/DMO/'
nout=80
info = load.info.Info(nout=nout, base=wdir, load=True)


# In[ ]:

s = load.sim.Sim(nout, wdir, setup=True)
s.add_part(ptypes=["dm id pos mass"], load=True)


# In[14]:

hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="RS", info=info)
hall.load()
region = smp.set_region(xc=0.5, yc=0.5, zc=0.5, radius=0.5) 


fig = plt.figure()
ax1 = fig.add_subplot(111)

#ptimg = img.ptden2d
#img.ptden2d.plot_2d_den(save=False, show=False, vmin=1e13, vmax=1e20, dpi=200, axes=ax1)

pp.pp_halo(hall, 800, region=region, axes=ax1, rscale=40, name=True)


# In[15]:

plt.show()


# In[ ]:



