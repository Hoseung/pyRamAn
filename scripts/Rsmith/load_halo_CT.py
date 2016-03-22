# -*- coding: utf-8 -*-
"""
Created on Sat Feb 14 14:49:22 2015

@author: hoseung

PURPOSE:
    halo mass evolution catalog for Jinsu's phase space diagram study.

OUTPUT:
    Consists of 10 columns, for each galaxy in each snapshot.
    #  ID   x   y   z   vx   vy   vz   Rvir   Mvir  Mass
    ((nnouts*ngal) x ncolumns) number of values.
    +
    array of snapshot redshifts at the end of file.

SAMPLING:
    First, halos that end up within 'rvir' * the cluster virial radius are chosen.
    Among them, only those with complete main prg tree are finally selected.
    Additionally, halo mass cut (Mcut) is available. 
"""

import numpy as np

def filter_halo_mass(data, Mcut=None):
    """ Returns indices of halos more massive tha Mcut"""
    m = np.array(data['m'][0])
    ind =np.where(m > Mcut)[0]
    print("# of halos:",len(ind))
    return ind


def n_most_massive(data, massive_count=1000):
    """ Returns indicies of top N most massive halos,
    massive_count = 1000 by default """
    m = np.array(data['m'][0])
    i = np.argsort(m)
    ind = i[:-1 - massive_count:-1]
    return ind


def filter_halo_pnum(data, Ncut=1000):
    """ Returns indicies of halos with more than Ncut particles"""
    npart = np.array(data['np'][0])
    ind =np.where(npart > Ncut)[0]
    print("# of halos:",len(ind))
    return ind


def extract_halos_within(halo, ind_center, scale=1.0):
    '''
    Returns halos within SCALE * Rvir of the central halo.

    def extract_halos_within(halos, ind_center, scale=1.0)
    halos : halo finder output (single snapshot)
    ind_center : index of central halo
    scale : multiplying factor to the Rvir of the central halo
    '''
    import numpy as np
    import utils.sampling as smp

    xc = halo.data['x'][i_center]
    yc = halo.data['y'][i_center]
    zc = halo.data['z'][i_center]
    rvir= halo.data['rvir'][i_center]

    xx = halo.data['x']
    yy = halo.data['y']
    zz = halo.data['z']
    m = np.array(halo.data['mvir'])

    dd = smp.distance_to([xc,yc,zc],[xx,yy,zz])

    Mcut = 1e11
    i_m = m > Mcut
    i_ok = np.logical_and(dd < (rvir * scale), i_m)

    return i_ok


#%%
''' Cluster 05101, cluster subhaloes (at the final snapshot)
'''

from tree import treemodule
from tree import halomodule
import tree.treeutils as tru
import numpy as np
trees = treemodule.CTree()
wdir = '/home/hoseung/Work/data/AGN2/'
trees.load(filename= wdir + 'rhalo/rockstar_halos/trees/tree_0_0_0.dat')

options = ['N most massive', '>1e13', 'nDM']
option = options[1]
n_massive = 500
fixed_position = True
Ncut = 120
#%%
nout_ini = 50
nout_fi = 132
nout_ini_hal = 10
rvir=3.0
# nout_halo = 122 == nout 10
# nout_halo = 0   == nout 132
#%%
nouts = range(nout_fi, nout_ini, -1)
Nnouts = len(nouts)

#%%
try:
    f = open(wdir + 'satellite_halos.txt', 'w')
    f_properties = open(wdir + 'satellite_halos_prop.txt', 'w')
except:
    print("No filename is given.\n Try write_halo_xyz(x,y,z,r,filename = fn)")

try:
    f = open(wdir + 'satellite_halos.txt', 'w')
    f_properties = open(wdir + 'satellite_halos_prop.txt', 'w')
except:
    print("No filename is given.\n Try write_halo_xyz(x,y,z,r,filename = fn)")

#%%
# in case of CT,  final = 81 & initial = 0  
tfin = trees.data[np.where(trees.data['nout'] == nout_fi - nout_ini - 1)]
tini = trees.data[np.where(trees.data['nout'] == 0)]

halo = halomodule.Halo(base=wdir, nout=nout_fi, halofinder="RS")
halo.load()
#%%

i_center = np.where(halo.data['np'] == max(halo.data['np']))
i_satellites = extract_halos_within(halo, i_center, scale=rvir)[0]
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))
      
# halos found inside the cluster and has tree back to nout_ini
halo_list = halo.data['id'][i_satellites]
#print(halo_list)

#%%
import utils.util
utils.util.reimport(tru)
good_trees=[]
n_good = 0
for hid in halo_list:
    prgs, inds = tru.get_main_prg(trees, haloid=hid, unique_id=False)
    if prgs is False:
        print("No tree for {}".format(hid))
        continue
    n_good +=1
#    print(len(prgs))
    good_trees.append(inds)

# again
inds_arr = np.zeros((n_good, 82), dtype=int)
n_good = 0
for hid in halo_list:
    prgs, inds = tru.get_main_prg(trees, haloid=hid, unique_id=False)
    if prgs is False:
        print("No tree for {}".format(hid))
        continue
    inds_arr[n_good] = inds
    n_good +=1
#    print(trees.data['mvir'][inds[0]] / halo.data['mvir'][hid])
    # trees and halo refer to the same halo, but quantities are slightly different.
    # .. why?? 
    
#%%
import load
f.write(" #      ID        x          y         z[Mpc]       vx      vy     vz[km/s]")
f.write("    Rvir(Mpc)      Mvir      Mass[Msun]\n")
zred=[]
for inout, nout in enumerate(nouts):
    info = load.info.Info(nout=nout, base=wdir)
    info.read_info()
    ind = inds_arr[:, inout]
    x = trees.data['x'][ind]# * info.pboxsize
    y = trees.data['y'][ind]# * info.pboxsize
    z = trees.data['z'][ind]# * info.pboxsize
    vx = trees.data['vx'][ind]# * info.kms
    vy = trees.data['vy'][ind]# * info.kms
    vz = trees.data['vz'][ind]# * info.kms
    r = trees.data['rvir'][ind]# * info.pboxsize
    m = trees.data['mvir'][ind]# * 1e11 
    m2 = trees.data['m200b'][ind]
    ids = [int(i) for i in trees.data['Orig_halo_id'][ind]]

    for i in range(len(ids)):
        f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(i,ids[i],x[i],y[i],z[i]))
        f.write("  {:.3f}  {:.3f}  {:.3f}".format(vx[i],vy[i],vz[i]))
        f.write("  {:.6f}  {:.0f}  {:.0f}     \n".format(r[i],m[i],m2[i]))
    zred.append(str(info.zred))

f.write(" Redshifts  \n")
for i, nout in enumerate(nouts):
    f.write("{0} ".format(zred[i]))

f.close()