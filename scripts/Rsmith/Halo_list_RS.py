
# coding: utf-8

# In[1]:


"""
Created on Tue Mar  3 13:12:43 2015
write a text file for Blender visualization

x  y  z  r

x  y  z  r

x  y  z  r

In this script top N massive halos
@author: hoseung
"""

import numpy as np

def distance_to(xc, xx):
    import numpy as np
    return np.sqrt([(xc[0] - xx[0])**2 + (xc[1] - xx[1])**2 + (xc[2] - xx[2])**2])


def extract_halos_within(halos, i_center, scale=1.0):
    import numpy as np
    '''
    Returns halos within SCALE * Rvir of the central halo.

    def extract_halos_within(halos, ind_center, scale=1.0)
    halos : halo finder output (single snapshot)
    ind_center : index of central halo
    scale : multiplying factor to the Rvir of the central halo
    '''
    xc = halos['p'][0][0][i_center]
    yc = halos['p'][0][1][i_center]
    zc = halos['p'][0][2][i_center]
    rvir= halos['rvir'][0][i_center]

    xx = halos['p'][0][0]
    yy = halos['p'][0][1]
    zz = halos['p'][0][2]
    m = np.array(data['m'][0])

    dd = distance_to([xc,yc,zc],[xx,yy,zz])

    Mcut = 1e11
    i_m = m > Mcut
    i_ok = np.logical_and(dd < (rvir * scale), i_m)

    return i_ok


''' Cluster 05101, cluster subhaloes (at the final snapshot)
'''
n_massive = 500
include_id = False
fixed_position = True
Ncut = 120
work_dir = '/home/hoseung/Work/data/AGN2/'
nout_ini = 131
nout_fi = 132
nouts = range(nout_fi, nout_ini, -1)
Nnouts = len(nouts)


try:
    f = open(work_dir + 'satellite_halo_trees.txt', 'w')
except:
    print("No filename is given.\n Try write_halo_xyz(x,y,z,r,filename = fn)")

import tree.treeutils as tru
import tree.load_RShalo as rsh
# get_main_tree
# gal_list
import pickle

dir_halo = work_dir + "rhalo/rockstar_halos/"
f_halo = dir_halo + "tree_py.dat"

# Open the file and call pickle.load.
with open(f_halo, "rb") as f_halo:
    data = pickle.load(f_halo)


# Gals = satellite halos inside the zoomed-in cluster above a mass cut.
all_final_halo = tru.gal_list(data)
# tru.gal_list returns the list of galaxies at the final snapshot.
i_center = tru.get_center(data)
# No # particle information. Devise a new way.

fn_halo = dir_halo + 'halos_py/halos_82.py'
## There is no halos_82.py.
#  You have to convert ascii files into .py or .pickle before.

#tree_final = data[i_final]
tree_final = rsh.read_halo_all(fn_halo)
##

i_satellites = extract_halos_within(data, i_center, scale=1.0)[0]
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))


cnt=150
ngal = len(gals)
#cnt = range(ngal)

# loop over individual galaxies
for thisgal in gals[cnt:cnt+1]:
    print(thisgal)
    if (cnt % 10 == 0): print(cnt)
    tree = tru.get_main_tree(data, thisgal)
    ind_tree = np.zeros(len(tree), dtype=np.int)
    for i in range(sum(x > 0 for x in tree)): ind_tree[i] = np.where(data['id'] == tree[i])[0]
    # Why first element?
    tree = data[ind_tree] # td = tree

    x = data['p'][0][0][i_satellites]
    y = data['p'][0][1][i_satellites]
    z = data['p'][0][2][i_satellites]
    r = data['rvir'][0][i_satellites]

    if include_id is True:
         dd = np.column_stack([x, y, z, r, data['hnu'][0][ind]])
         for i in range(dd.shape[0]):
             f.write("{0}   {1}   {2}   {3}  {4} \n".format(
                     dd[i][0],dd[i][1],dd[i][2],dd[i][3],int(dd[i][4])))
    else:
        dd = np.column_stack([x, y, z, r])
        for i in range(dd.shape[0]):
             f.write("{0}   {1}   {2}   {3} \n".format(
                     dd[i][0],dd[i][1],dd[i][2],dd[i][3]))

f.close()

