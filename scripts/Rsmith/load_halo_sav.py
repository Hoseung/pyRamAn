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


def extract_halos_within(halos, ind_center, scale=1.0, Mcut=1e4):
    '''
    Returns halos within SCALE * Rvir of the central halo.

    def extract_halos_within(halos, ind_center, scale=1.0)
    halos : halo finder output (single snapshot)
    ind_center : index of central halo
    scale : multiplying factor to the Rvir of the central halo
    '''
    import numpy as np
    import utils.sampling as smp

    xc = halos['p'][0][0][i_center]
    yc = halos['p'][0][1][i_center]
    zc = halos['p'][0][2][i_center]
    rvir= halos['rvir'][0][i_center]

    xx = halos['p'][0][0]
    yy = halos['p'][0][1]
    zz = halos['p'][0][2]
    m = np.array(halos['m'][0])

    dd = smp.distance_to([xc,yc,zc],[xx,yy,zz])

    i_m = m > Mcut
    i_ok = np.logical_and(dd < (rvir * scale), i_m)

    return i_ok


#%%
''' Cluster 05101, cluster subhaloes (at the final snapshot)
'''

import tree.hmutil as hmu

options = ['N most massive', '>1e13', 'nDM']
option = options[1]
n_massive = 500
fixed_position = True
Ncut = 120
Mcut = 1e7
#work_dir = '/home/hoseung/Work/data/AGN2/'
work_dir = './'
nout_ini = 30
nout_fi = 132
nout_ini_hal = 10
rvir=3.0
# nout_halo = 122 == nout 10
# nout_halo = 0   == nout 132
nouts = range(nout_fi, nout_ini -1, -1)
Nnouts = len(nouts)

try:
    f = open(work_dir + 'satellite_halos.txt', 'w')
#    f_properties = open(work_dir + 'satellite_halos_prop.txt', 'w')
except:
    print("No filename is given.\n Try write_halo_xyz(x,y,z,r,filename = fn)")

from tree import tmtree
tree = tmtree.load(work_dir=work_dir, filename="halo/TMtree.fits")
tfin = tree[np.where(tree['NOUT'] == 0)]
tini = tree[np.where(tree['NOUT'] == nout_fi - nout_ini)]

halo = hmo.Halo(nout=nout_fi, base=work_dir, halofinder="HM", load=True)
#%%
i_center = np.argmax(halo.data['np'])
import utils.sampling as smp
import tree.halomodule as hmo
i_satellites = smp.extract_halos_within(halo.data, i_center, scale=rvir, Mcut=Mcut)
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))

# halos found inside the cluster and has tree back to nout_ini
halo_list = halo.data['id'][i_satellites]
#print(halo_list)
h_ind_ok, halo_ok = tmtree.check_tree_complete(tree, 0, nout_fi - nout_ini, halo_list)
print(len(halo_ok))

from utils import match
import load
f.write(" #      ID        x          y         z[Mpc]       vx      vy     vz[km/s]")
f.write("    Rvir(Mpc)      Mvir      Mass[Msun]\n")
zred=[]
for inout, nout in enumerate(nouts):
    info = load.info.Info()
    info.setup(nout=nout, base=work_dir)
    info.read_info()
#    data = hmu.load_data(nout, work_dir=work_dir, normalize=True) # load .sav halo file and normalize it to code unit.
    halo = hmo.Halo(nout=nout_fi, base=work_dir, halofinder="HM", load=True, info=info)
#    fname = work_dir + 'halos_py/halos_' + '031' + '.pickle'
#    data = load_halo_py(fname)
    ind = match.match_list_ind(data['id'], halo_ok[:,inout])

    x = halo.data['x'][ind] * info.pboxsize
    y = halo.data['y'][ind] * info.pboxsize
    z = halo.data['z'][ind] * info.pboxsize
    vx = halo.data['vx'][ind]# * info.kms
    vy = halo.data['vy'][ind]# * info.kms
    vz = halo.data['vz'][ind]# * info.kms
    r = halo.data['rvir'][ind] * info.pboxsize
    m = halo.data['mvir'][ind] * 1e11 
    m2 = halo.data['m'][ind]
    ids = [int(i) for i in halo.data['id'][ind]]

    for i in range(len(ids)):
        f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(i,ids[i],x[i],y[i],z[i]))
        f.write("  {:.3f}  {:.3f}  {:.3f}".format(vx[i],vy[i],vz[i]))
        f.write("  {:.6f}  {:.0f}  {:.0f}     \n".format(r[i],m[i],m2[i]))
    zred.append(str(info.zred))

f.write(" Redshifts  \n")
for i, nout in enumerate(nouts):
    f.write("{0} ".format(zred[i]))

f.close()