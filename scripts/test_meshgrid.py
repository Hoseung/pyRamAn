# -*- coding: utf-8 -*-
"""
Created on Thu May 28 16:19:52 2015

Mesh grid test
@author: hoseung
"""

def get_sigma(vel):
    import numpy as np
    return np.std(vel)
    

# Load data
import h5py as hdf
wdir = '/home/hoseung/Work/data/AGN2/'
infile = hdf.File(wdir+'005282gal.hdf5', "r")

x = infile['galaxy/star/x'].value
y = infile['galaxy/star/y'].value # There must be a way to access several at once!
z = infile['galaxy/star/z'].value

m = infile['galaxy/star/m'].value




ar/vx'].value
vy = infile['galaxy/star/vy'].value # There must be a way to access several at once!
vz = infile['galaxy/star/vz'].value

#%%
import matplotlib.pyplot as plt
import numpy as np

# NGP
npix = 50
nx = npix
ny = npix

xx = (x - min(x)) / x.ptp() * nx
yy = (y - min(y)) / y.ptp() * ny

ngx = np.clip(np.fix(xx), 0, nx-1)
ngy = np.clip(np.fix(yy), 0, ny-1)

indices = ngx + ngy * nx

# 
sigmap=np.zeros(nx * ny, dtype=float)
for i in range(nx * ny):
    ind = np.where(indices == i)
    sigmap[i]=get_sigma(vz[ind])


# projected quantity
vmap = np.zeros(nx * ny, dtype=float)
for i, ind in enumerate(indices):
#    fimg[ind] = vz[i]
#    fimg[ind] = m[i]
    vmap[ind] += m[i]

plt.imshow(vmap.reshape((nx, ny)).squeeze())
cbar = plt.colorbar()
#cbar.set_label(r'$Log\Sigma$ $(M_{\odot} kpc^{-2})$',  rotation=270, labelpad=10)
plt.savefig(wdir + 'vmap.png')


plt.imshow(sigmap.reshape((nx, ny)).squeeze())
cbar = plt.colorbar()
plt.savefig(wdir + 'sigmap.png')
#%%