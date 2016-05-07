# -*- coding: utf-8 -*-
"""
Created on Sat Feb 14 14:49:22 2015

@author: hoseung

Complete. 
Directly read treebricks file. 

PURPOSE:
    write a text file for Jinsu's phase diagram study.

OUTPUT:
    Consist of multiple columns for each galaxy in each snapshot.
    #  ID   x   y   z   vx   vy   vz   Rvir   Mvir  Mass
    ((nnouts*ngal) x ncolumns) number of values.
    +
    array of snapshot redshifts at the end of file.

SAMPLING:
    First, halos that end up in the main cluster are chosen.
    Among them, only those with complete main prg tree are finally selected.


"""

#from scipy.io.idl import readsav
import numpy as np

def get_idx(tree, hnus, nout=None):
    i_nout = np.where(tree.field('NOUT') == nout)
    i_halo = match_list_ind(tree[i_nout].field('HALNUM'), hnus)

    return tree[i_nout[i_halo]].field('IDX')


def filter_halo_mass(data, Mcut=None):
    m = np.array(data['m'][0])
    #ind = m > Mcut
    #print("# of halos:",sum(ind))
    ind =np.where(m > Mcut)[0]
    print("# of halos:",len(ind))
    return ind


def n_most_massive(data, mass_count=1000):
    m = np.array(data['m'])
    i = np.argsort(m)
    ind = i[:-1 - mass_count:-1]
    return ind


def filter_halo_pnum(data, Ncut=1000):
    npart = np.array(data['np'])
    ind =np.where(npart > Ncut)[0]
    print("# of halos:",len(ind))
    return ind


def extract_halos_within(halos, ind_center, scale=1.0, Mcut=1e5):
    import numpy as np
    import utils.sampling as smp
    '''
    Returns indices of halos within SCALE * Rvir of the central halo.

    def extract_halos_within(halos, ind_center, scale=1.0)
    halos : halo finder output (single snapshot)
    ind_center : index of central halo
    scale : multiplying factor to the Rvir of the central halo
    '''
    xc = halos['x'][i_center]
    yc = halos['y'][i_center]
    zc = halos['z'][i_center]
    rvir= halos['rvir'][i_center]

    xx = halos['x']
    yy = halos['y']
    zz = halos['z']
    m = np.array(halos['m'])

    dd = smp.distance_to([xc,yc,zc],[xx,yy,zz])

    return (dd < rvir * scale) * (m > Mcut)


#%%
''' Cluster 05101, cluster subhaloes (at the final snapshot)
'''

import tree.halomodule as hmo 
from utils import match
import load

fixed_position = True
#work_dir = './'
work_dir = '/home/hoseung/Work/data/05427/'
nout_ini = 37
# 27 : z=4
# 37 : z=3
# 20 : ~ z=5
nout_fi = 187
rvir=3.0
# nout_halo = 122 == nout 10
# nout_halo = 0   == nout 132
nouts = range(nout_fi, nout_ini -1, -1) 
Nnouts = len(nouts)

from tree import tmtree
tree = tmtree.load(work_dir=work_dir, filename="halo/TMtree.fits")
tfin = tree[np.where(tree['NOUT'] == 0)]
tini = tree[np.where(tree['NOUT'] == nout_fi - nout_ini)]
#%%
info = load.info.Info()
info.setup(nout=nout_fi, base=work_dir)
info.load()
hh = hmo.Halo(base=work_dir, nout=nout_fi, halofinder='HM', info=info)
hh.load()

#halo = hmu.load_data(nout_fi, work_dir=work_dir, normalize=True)
i_center = np.where(hh.data['np'] == max(hh.data['np']))[0]
i_satellites = extract_halos_within(hh.data, i_center, scale=3.0)
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))

# halos found inside the cluster and has tree back to nout_ini
halo_list = hh.data['id'][i_satellites]
#print(halo_list)
h_ind_ok, halo_ok = tmtree.check_tree_complete(tree, 0, nout_fi - nout_ini, halo_list)
print(len(halo_ok))


try:
    f = open(work_dir + 'satellite_halos.txt', 'w')
except:
    print("No filename is given.\n Try write_halo_xyz(x,y,z,r,filename = fn)")

f.write(" #      ID        x          y         z[Mpc]       vx      vy     vz[km/s]")
f.write("    Rvir(Mpc)      Mvir      Mass[10e10Msun]\n")
zred=[]
# It works, but too slow!
for inout, nout in enumerate(nouts):
    print(nout)
    info = load.info.Info()
    info.setup(nout=nout, base=work_dir)
    info.load()

    fn = work_dir + 'halo/tree_bricks' + str(nout).zfill(3)
    hh = hmo.Halo(base=work_dir, nout=nout, halofinder='HM', info=info)
    hh.load()
#    hh = tree.hmhalo.Hmhalo(fn)
    data = hh.data

    ind = match.match_list_ind(data['id'], halo_ok[:,inout])

    x = data['x'][ind] * info.pboxsize
    y = data['y'][ind] * info.pboxsize
    z = data['z'][ind] * info.pboxsize
    vx = data['vx'][ind]# * info.kms
    vy = data['vy'][ind]# * info.kms
    vz = data['vz'][ind]# * info.kms
    r = data['rvir'][ind] * info.pboxsize
    m = data['mvir'][ind]
    m2 = data['m'][ind]
    ids = [int(i) for i in data['id'][ind]]

    for i in range(len(ids)):
        f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(i,ids[i],x[i],y[i],z[i]))
        f.write(" {1:>.6f}  {:.3f}  {:.3f}".format(vx[i],vy[i],vz[i]))
        f.write("  {:.6f}  {:.0f}  {:.0f}     \n".format(r[i],m[i],m2[i]))
    zred.append(str(info.zred))

f.write(" Redshifts  \n")
for i, nout in enumerate(nouts):
    f.write("{0} ".format(zred[i]))

f.close()
