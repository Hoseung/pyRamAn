# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 22:15:31 2015

This routine makes plots to test CT outputs.
Test halo porperty consistency.

@author: hoseung
"""
import numpy as np
import tree.treeutils as tru
import tree.treeplots as trp
import tree.rshalo as rsh
import utils.sampling as smp
from utils import util
import pickle


# Load tree
wdir = '/home/hoseung/Work/data/05427/'
dir_halo = wdir + "rhalo/rockstar_halos/"

fn_halo = dir_halo + 'halos_81.ascii'

nout_ini = 0
nout_fi = 81
nouts = range(nout_fi, nout_ini, -1)
Nnouts = len(nouts)

# RS tree
f_tree = wdir + "rhalo/Trees/tree.pickle"
with open(f_tree, "rb") as ft:
    rstree = pickle.load(ft)

from astropy.io import fits
from astropy.table import Table
import tree
import load

info = load.info.Info(nout=132, base=wdir)

# HM tree
data = fits.getdata(wdir + "halo/TMtree.fits", 1)
hmtree = Table(data)
#%%

# Gals = satellite halos inside the zoomed-in cluster above a mass cut.
# tru.gal_list returns the list of galaxies at the final snapshot.
# Pick some of trees (or all)
#final_halos_rs = tru.final_halo_list(rstree)
i = np.where(rstree['nout'] == 81) # 1 ~ 82,...-> 0 ~ 81
x = rstree['x'][i]  * info.pboxsize / 200
y = rstree['y'][i]  * info.pboxsize / 200
z = rstree['z'][i]  * info.pboxsize / 200

#%%
#final_halos_hm = tru.final_halo_list(hmtree)
inout = np.where(hmtree['NOUT'] == 0)[0]
final_halos_hm = hmtree["HALNUM"][inout]   # 8150
#print(len(inout))
x2 = hmtree["XP"][inout,0] + info.pboxsize * 0.5  # in Mpc
y2 = hmtree["XP"][inout,1] + info.pboxsize * 0.5
z2 = hmtree["XP"][inout,2] + info.pboxsize * 0.5


#%%
import matplotlib.pylab as plt

# find matching halos and build common halo list (2-D)
# x is larger array
#def match_diff_arr(x1, x2, y1, y2, z1=None, z2=None, tolerance=None):
# x[i] = x2[mi[i]]
# if no match, mi[i] = -1
import zoom.utils as zut

match_ind = zut.match_diff_arr(y2, y, x2, x, z1 = z2, z2 = z, tolerance=0.002) # 100kpc
# How can I add mass?

#inda = np.where(match_ind >= 0)[0]
indb = match_ind[match_ind >= 0]


#%%
#x = x[inda]
#x2 = x2[indb]
plt.close()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(z[indb])
ax.plot(z2[inda]+0.5)
#ax.hist(z)
#ax.hist(z2, range=(110, 140))
#%%
for i,j in zip(x[0:10], x2[0:10]):
    print(i,j)

#%%
for i in range(sum(match_ind >= 0)):
    common_halos[0]

## ID list

#%%
# mass-produce plots of halo properties.
quantities=["sam_mvir","mvir","rvir","rs","vrms","vmax"
    ,"jx","jy","jz","spin","m200b","m200c","m500c","m2500c"
    ,"xoff","voff","btoc","ctoa","ax","ay","az"]

normalizer=np.array([1e-11,1e-11,1,1,1,1,1,1,1e-11
        ,1,1e-11,1e-11,1e-11,1e-11,1,1,1,1,1,1,1])

for i, hal in enumerate(zip(final_halos_rs, final_halos_hm)):
    print(i, hal)
    thistree = tru.get_main_prg_tree(rstree, hal)

    fn_save = str(hal) + 'halo_all.pdf'
#    trp.plot_all(tree, hal, save=True, out_dir=work_dir, fn_save=fn_save,
#                 nrows=4, ncols=math.ceil(len(quantities)/4),
#                 quantities=quantities, normalizer=normalizer)
#    trp.plot_all_multiPDF(tree, hal, out_dir=work_dir + 'RS_trees/', fn_save=fn_save,
#                 nrows=2, ncols=2,
#                 quantities=quantities, normalizer=normalizer)
    trp.trajectory3D(thistree, hal, save=wdir + 'RS_trees/')
# # add position plots in code unit (comoving scale)

# One ps per halo.


# 4 figures per page.