# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 12:43:15 2015

Halo trace + halo Evolution history

@author: hoseung
"""

#%%
import tree.hmutil as hmu
import tree.halomodule as hmo
import numpy as np
import utils.util
from tree import tmtree
import utils.sampling as smp

options = ['N most massive', '>1e13', 'nDM']
option = options[1]
n_massive = 500
fixed_position = True
Ncut = 120
Mcut = 1e7
work_dir = '/home/hoseung/Work/data/05427/'
#work_dir = './'
nout_ini = 37
nout_fi = 187
nout_ini_hal = 10
rvir=3.0
# nout_halo = 122 == nout 10
# nout_halo = 0   == nout 132
nouts = range(nout_fi, nout_ini -1, -1)
Nnouts = len(nouts)

tree = tmtree.load(work_dir=work_dir, filename="halo/TMtree.fits")
tfin = tree[np.where(tree['NOUT'] == 0)]
tini = tree[np.where(tree['NOUT'] == nout_fi - nout_ini)]

utils.util.reimport(hmo)
#halo = hmu.load_data(nout_fi, work_dir=work_dir, normalize=True)
hh = hmo.Halo(nout=nout_fi, base=work_dir, halofinder="HM", load=True)

#%%
import draw.pp
import matplotlib.pyplot as plt

i_center=1599
hh = hmo.Halo(nout=nout_fi, base=work_dir, halofinder="HM", load=True)
region = smp.set_region(xc = hh.data['x'][i_center],
                        yc = hh.data['y'][i_center],
                        zc = hh.data['z'][i_center],
                        radius = 0.001)
#                        radius = hh.data['rvir'][i_center] * 5)
ind_prg, a_prg = tmtree.get_main_prg(tree, i_center, nout_ini=nout_fi - nout_ini)
plt.ioff()
#%%
npix = 400
for i, nout in enumerate(range(nout_fi,nout_ini,-1)):
    hh = hmo.Halo(nout=nout, base=work_dir, halofinder="HM", load=True)
    i_hal_now = np.where(hh.data['id'] == a_prg[0,i])[0]
    i_satellites = smp.extract_halos_within(hh.data,
                                    i_hal_now,
                                    scale=0.001/hh.data['rvir'][i_hal_now],
                                    Mcut=1e10)

    region = smp.set_region(xc = hh.data['x'][i_hal_now],
                        yc = hh.data['y'][i_hal_now],
                        zc = hh.data['z'][i_hal_now],
                        radius = 0.001)

    fig, ax = plt.subplots()
    draw.pp.pp_halo(hh, npix,
                    ind=np.arange(len(i_satellites))[i_satellites],
                    axes=ax,
                    name=True,
                    region=region)

    ax.set_xlim([-0.5 * npix, 1.5 * npix])
    ax.set_ylim([-0.5 * npix, 1.5 * npix])
    ax.set_aspect('equal', 'datalim')
    #plt.show()
    plt.savefig(work_dir + 'galaxy_' + str(i_center).zfill(5) + '/' + str(nout) + '.png', dpi=200)
    plt.close()

#%%
npix = 400
# Version 1
# All halos move w.r.t. the final halo position. 
region = smp.set_region(xc = hh.data['x'][i_center],
                        yc = hh.data['y'][i_center],
                        zc = hh.data['z'][i_center],
                        radius = 0.001)

for i, nout in enumerate(range(nout_fi,nout_ini,-1)):
    hh = hmo.Halo(nout=nout, base=work_dir, halofinder="HM", load=True)
    i_hal_now = np.where(hh.data['id'] == a_prg[0,i])[0]
    i_satellites = smp.extract_halos_within(hh.data,
                                    i_hal_now,
                                    scale=0.001/hh.data['rvir'][i_hal_now],
                                    Mcut=1e10)
    fig, ax = plt.subplots()
    draw.pp.pp_halo(hh, npix,
                    ind=np.arange(len(i_satellites))[i_satellites],
                    axes=ax,
                    name=True,
                    region=region)

    ax.set_xlim([-1 * npix, 6 * npix])
    ax.set_ylim([-1 * npix, 6 * npix])
    ax.set_aspect('equal', 'datalim')
    #plt.show()
    plt.savefig(work_dir + str(nout) + '.png', dpi=200)
    plt.close()
#ax.set_autoscale_on(False)



#anim = animation.FuncAnimation(fig, animate, frames=30)
#anim.save('demoanimation.gif', writer='imagemagick', fps=4);



#%%
try:
    f = open(work_dir + 'satellite_halos.txt', 'w')
#    f_properties = open(work_dir + 'satellite_halos_prop.txt', 'w')
except:
    print("No filename is given.\n Try write_halo_xyz(x,y,z,r,filename = fn)")


#%%
#i_center = np.argmax(halo.data['np'])
i_center=1591


i_satellites = smp.extract_halos_within(halo.data, i_center, scale=rvir, Mcut=Mcut)
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))

# halos found inside the cluster and has tree back to nout_ini
halo_list = halo.data['id'][i_satellites]
#print(halo_list)
#h_ind_ok, halo_ok = tmtree.check_tree_complete(tree, 0, nout_fi - nout_ini, halo_list)

# indicies of all nearby halos to the target halo.
ind_prg, a_prg = tmtree.get_main_prg(tree, halo_list)
# Tree nout goes opposite. Need to be corrected. 


