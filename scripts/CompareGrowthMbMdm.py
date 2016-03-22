# -*- coding: utf-8 -*-
"""
Created on Sun Feb 15 17:15:05 2015

Compare growths of baryon mass and DM mass of galaxies to test if it is valid to
assume that the galaxy baryon component grow propotional to the DM halo mass.
(SAM!)

Use Rockstar halo outputs!!

@author: hoseung
"""
work_dir="/home/hoseung/Work/data/AGN2/"

gal_dir = work_dir + "galaxies132/galaxies"

halo_dir = work_dir + "halos_py/"

#fn_tree = work_dir + "rhalo/rockstar_halos/trees/tree_0_0_0.dat"
fn_tree = work_dir + "rhalo/rockstar_halos/tree_py.dat"


# get galaxy evolution tree
import tree
import pickle

#%%
#  tt = pickle.load(halo_dir + ) It's tree..
#  tree file is in .dat form. one in original CT output,
#  Another is in python-saved form.
f_pickle = open(fn_tree, 'rb') # Open file first and pass to pickle. .. WHY??
tt = pickle.load(f_pickle)


#%%



# read galaxy at each nout
# Use HaloMaker output for the moment.
# Which in turn requires 1:1 match of HaloMaker halos and Rockstar halos
# (using a separate script)




# Measure mass of young star
# Measure stellar mass change
# Measure halo mass change

# build an array of Mstellar_arr and Mhalo_arr

# plot