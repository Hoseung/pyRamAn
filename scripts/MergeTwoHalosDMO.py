# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 21:15:33 2015

Load two different halos from DMO runs and merge them.

@author: hoseung
"""

from tree import halomodule
from utils import util
import utils.match as mtc
import load
import numpy as np
import matplotlib.pyplot as plt
import zoom

nout = 80

work_dir = '/home/hoseung/Work/data/DMO/'

info = load.info.Info(nout=nout, base=work_dir)
info.read_info()
h = halomodule.Halo(base=work_dir, nout=nout, info=info)

h.set_halofinder('HaloMaker')
h.load_hm()

# Merger two simulation with offset.
work_dir2 = '/home/hoseung/Work/data/DMO_000/'

info2 = load.info.Info(nout=nout, base=work_dir2)
info2.read_info()
h2 = halomodule.Halo(base=work_dir2, nout=nout, info=info2)
h2.set_halofinder('HaloMaker')
h2.load_hm()


#%%
mass_cut = 6e13
scale = 10

xr = np.array([110000, 190000]) / info.cboxsize / 1000.
yr = np.array([110000, 190000]) / info.cboxsize / 1000.
zr = np.array([110000, 190000]) / info.cboxsize / 1000.


ind = np.where((h.data.mvir > mass_cut) &
                (h.data.x > xr[0]) & (h.data.x < xr[1]) &
                (h.data.y > yr[0]) & (h.data.y < yr[1]) &
                (h.data.z > zr[0]) & (h.data.z < zr[1]) )

x = h.data.x[ind] * info.cboxsize * 1000.
y = h.data.y[ind] * info.cboxsize * 1000.
z = h.data.z[ind] * info.cboxsize * 1000.
r = h.data.rvir[ind] * 0.5 * info.cboxsize * 1000. * scale

# need a reference to an axes to keep plotting on the same plot.
fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.scatter(x,y,r, facecolor='none', edgecolor='b', label='555')
# size of scatter plot marker scales with the figure scale.
# Convinient but unphysical.

from draw import pp
pp.circle_scatter(ax1, x, y, r, facecolor='none', edgecolor='b', label='555')

#%%
xr2 = np.array([10000, 90000]) / info.cboxsize / 1000.
yr2 = np.array([10000, 90000]) / info.cboxsize / 1000.
zr2 = np.array([10000, 90000]) / info.cboxsize / 1000.


ind2 = np.where((h2.data.mvir > mass_cut) &
                (h2.data.x > xr2[0]) & (h2.data.x < xr2[1]) &
                (h2.data.y > yr2[0]) & (h2.data.y < yr2[1]) &
                (h2.data.z > zr2[0]) & (h2.data.z < zr2[1]) )

x2 = h2.data.x[ind2] * info.cboxsize * 1000. + 100000
y2 = h2.data.y[ind2] * info.cboxsize * 1000. + 100000
z2 = h2.data.z[ind2] * info.cboxsize * 1000. + 100000
r2 = h2.data.rvir[ind2] * info.cboxsize * 1000. * scale

# square size = length of a side. but in circle_scatter it's the radius.
# So multiply by 2.!
pp.square_scatter(ax1, x2, y2, r2, facecolor='none', edgecolor='r', label='000')
#ax1.scatter(x2,y2,r2, facecolor='none', edgecolor='r', label='000', marker='s')


ax1.set_aspect("equal")
ax1.set_xlim(100000, 200000)
ax1.set_ylim(100000, 200000)

for i in range(len(x2)):
    ax1.annotate(str(i), (x2[i],y2[i]))
#plt.legend(loc='upper left')
plt.show()

#%%
# Matched indices 
mi = zoom.utils.match_diff_arr(x, x2, y, y2, tolerance=500)

for i in range(len(x)):
    print(zoom.utils.distance3d(x[i], y[i], z[i], x2[mi[i]], y2[mi[i]], z2[mi[i]]))
    