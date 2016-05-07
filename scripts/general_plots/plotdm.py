# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 01:08:05 2015

@author: hoseung
"""

#%%
import numpy as np
import load
import time
import utils.util as utl
import utils.prettyplot as ptt

base = '/home/hoseung/Work/data/group3826/'
nout=135
#base = '/home/hoseung/Work/data/Lmax_136/'
#nout=136
#base = '/home/hoseung/Work/data/C04466/'
#nout=181
#base = '/home/hoseung/Work/data/AGN2/'
#nout=51

foutput = 'field_Lmax136.npy'

import utils.sampling as smp

#region = smp.set_region(xc=0.205, yc=0.78, zc=0.395, radius=0.01)
#print(region)
#xrs = region["ranges"]
#xrs = [[0.21, 0.22], [0.77, 0.78], [0.39, 0.40]]
s = load.sim.Sim(nout, base) # , ranges=xrs)
#s.set_ranges(xrs)
#%%
s.add_part('dm')
print(" add_part ")
#s.part.dm['px'].ptp()

#region = smp.set_region(centers=[0.6,0.4,0.5],radius=0.005)
#print(" set_region ")

s.part.load()
#print(" set_part ")
#%%
region = s.part.search_zoomin(scale=0.5, load=True)


s.set_ranges(region["ranges"])
# sim.set_ranges does not automatically change sim.part.ranges
# Can I make sim.part.ranges to be t

print(" ------cc ")
s.part.load()
#%%

# part
x = s.part.dm['px']
y = s.part.dm['py']  # These are views. right?
z = s.part.dm['py']
m = s.part.dm['m'] * s.info.msun # part must be normalized already!

npix = 100

from draw import pp
import matplotlib.pyplot as plt
field = pp.den2d(x, y, z, m, npix, s.info, cic=True, norm_integer=True)
np.save(foutput, field)
# Save



#%%
import matplotlib.pyplot as plt
field = np.random.rand(300,300)
from matplotlib.colors import LogNorm
cname = 'brg'
vmin = field.max() * 1e-6
vmax = field.max()

p = plt.imshow(field, cmap=plt.get_cmap(cname), norm=LogNorm(vmin=vmin, vmax=vmax))
"""
nticks = 2
xr = np.asarray(region['xr'])
yr = np.asarray(region['yr'])
x_old, x_new = ptt.tickInPboxInt(xr * s.info.pboxsize, field.shape[0], nticks=nticks)
print(xr * s.info.pboxsize)
print(x_new)
print(x_old * npix)
y_old, y_new = ptt.tickInPboxInt(yr * s.info.pboxsize, field.shape[1], nticks=nticks)

p.xticks(x_old * npix, x_new)
p.yticks(y_old * npix, y_new)
"""

#legend
cbar = plt.colorbar(p)
cbar.ax.set_yticklabels(['0','1','2','>3'])
cbar.set_label(r'$Log\Sigma$ $(M_{\odot} kpc^{-2})$', rotation=270, labelpad=25)
plt.show()
plt.savefig('test.png', dpi=72)
plt.savefig('test.eps', dpi=200)
#%%
plt.close()

#plt.hist(field.flatten(), 50, normed=1, facecolor='green', alpha=0.75)
#plt.show()
# flatten 2D array before passing to .hist()

