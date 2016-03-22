# -*- coding: utf-8 -*-
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
import tree.treeutils as tru
import utils.sampling as smp
from utils import match
import pickle

''' Cluster 05101, cluster subhaloes (at the final snapshot)
'''

work_dir = '/home/hoseung/Work/data/AGN2/'
dir_halo = work_dir + "rhalo/rockstar_halos/"
f_tree = work_dir + "rhalo/tree.pickle"
#fn_halo = dir_halo + 'halos_81.ascii'
fn_halo = work_dir + 'halos_py/halos_081.pickle'
f_list = work_dir + 'satellite_halo_trees.txt'
f_property = work_dir + 'satellite_halo_trees_p.txt'
include_properties = False


nout_ini = 0
nout_fi = 81
nouts = range(nout_fi, nout_ini, -1)
Nnouts = len(nouts)


try:
    f = open(f_list, 'w')
    f_p = open(f_property, 'w')
except:
    print("No filename is given.\n Try write_halo_xyz(x,y,z,r,filename = fn)")


# Open the file and call pickle.load.
with open(f_tree, "rb") as ft:
    all_tree = pickle.load(ft)

# Gals = satellite halos inside the zoomed-in cluster above a mass cut.
# tru.gal_list returns the list of galaxies at the final snapshot.
all_final_halo = tru.final_halo_list(all_tree)

#halo_final = rsh.read_halo_all(fn_halo) # read ascii halo output
with open(fn_halo, 'rb') as f:
    halo_final = pickle.load(f)

#%%
i_center = np.where(halo_final['num_p'] == halo_final['num_p'].max())

i_satellites = smp.extract_halos_within(halo_final, i_center, scale=1.0)
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))

#%%
# halos found inside the cluster with complete tree
halo_list = halo_final['id'][i_satellites]
h_ind_ok, halo_ok = tru.check_tree_complete(all_tree, halo_list,
                                            nout_ini=nout_ini, nout_fi=nout_fi
                                            , return_idx=True)
print("{0} cluster halos have complete trees".format(len(halo_ok)))
#%%

# loop over individual galaxies
for inout, nout in enumerate(nouts):
    f_halo = work_dir + 'halos_py/halos_' + str(nout).zfill(3) + '.pickle'
    with open(f_halo, "rb") as f_hal:
        data = pickle.load(f_hal)
    ind = match.match_list_ind(data['id'], halo_ok[:,inout])
    x = data['x'][ind]
    y = data['y'][ind]
    z = data['z'][ind]
    r = data['rvir'][ind]

    if include_properties is True:
         ddp = np.column_stack([r, data['id'][0][ind]])
         for i in range(ddp.shape[0]):
             f_p.write("{0}   {1}   {2}   {3}  {4} \n".format(
                     ddp[i][0],ddp[i][1],ddp[i][2],ddp[i][3],int(ddp[i][4])))
    dd = np.column_stack([x, y, z, r])
    for i in range(dd.shape[0]):
         f.write("a{0}   {1}   {2}   {3}  {4}\n".format(
                 i,dd[i][0],dd[i][1],dd[i][2],dd[i][3]))


f.close()
f_p.close()

