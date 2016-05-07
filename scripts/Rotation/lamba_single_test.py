# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 14:44:18 2015

@author: hoseung
"""


import matplotlib.pyplot as plt

# Load catalog
import pickle
import numpy as np

wdir = '/home/hoseung/Work/data/05427/'

def load_pickle(fname):
    with open(fname, 'rb') as f:
        return pickle.load(f)

cat = load_pickle(wdir + 'catalog/' + 'catalog187_1275.pickle')

#%%

# 3D trajectory
xx = cat['xc']
yy = cat['yc']
zz = cat['zc']


from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=plt.figaspect(2.))
ax = fig.add_subplot(2, 1, 1)
ax.plot(cat['rgal'])

ax = fig.add_subplot(2,1,2, projection='3d')
print(len(xx))
#for i in range(len(xx)):
ax.plot(xx,yy,zz)

plt.show()

#%%

from tree import tmtree
#import utils.sampling as smp
work_dir = '/home/hoseung/Work/data/05427/'
tt = tmtree.load(work_dir=work_dir, filename="halo/TMtree.fits")

nout = 187
gal = 1441
i = np.where((tt['NOUT'] == 187 - nout) * (tt["HALNUM"] == gal))

tree = tt["TREE"][i][0]
tree = tree[tree > 0]

gal_idx = tree[0]
#%%
from draw import pp
fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
def plot_halo(tt, idx, ax, npix=800, **kwargs):
    x, y, z = tt["XP"][idx] / tt["BOXSIZE"][idx] +   0.5
    r = np.asarray(tt["RVIR"][idx] / tt["BOXSIZE"][idx])
    pp.circle_scatter(ax, x, y, r, facecolor='none', **kwargs)


def plot_halo3d(tt, idx, ax, c, **kwargs):
    x, y, z = tt["XP"][idx] / tt["BOXSIZE"][idx] +   0.5
    r = np.asarray(tt["RVIR"][idx] / tt["BOXSIZE"][idx]) * 100000
    ax.scatter(x,y,z, s=r, c=c, edgecolors='none')
#ax.scatter(xs, ys, zs, c=c, marker=m)    

def get_next_branch(tt, idx, ax):
    from matplotlib import cm
#    i = tt["IDX"] == idx))
    tree = tt["TREE"][idx]#[0]
#    print(tree)
    tree= tree[tree > 0]
    for i, ii in enumerate(tree):
        edgc = cm.jet(tt["NOUT"][ii])
    
#        plot_halo(tt, ii, ax, edgecolor=edgc)
        plot_halo3d(tt, ii, ax, edgc)
        get_next_branch(tt, ii, ax)

    return True  
    
from mpl_toolkits.mplot3d import Axes3D    
ax = fig.add_subplot(111, projection='3d')
get_next_branch(tt, gal_idx, ax)
ax.set_xlim([0.495,0.515])
ax.set_ylim([0.3,0.32])

#%%    

#    return get_next_branch(tree)
for i in aa:
    if sum(i > 0) > 1:
        print(i[0:5])
#%%
for branch in tree:
    get_next_branch(tt, branch)

