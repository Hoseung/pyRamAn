# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 15:13:39 2015

@author: hoseung
"""

# Conver Consistent tree file to pickle file.
import tree.treeutils as tru
base = "/home/hoseung/Work/data/AGN2/"
dir_halo = base + "rockstar_halos/"

f_CT = dir_halo + "trees/tree_0_0_0.dat" # AGN2
data = tru.load(f_CT)

import pickle 
f = dir_halo + "tree_py.dat"

# Open the file and call pickle.dump.
with open("f.pickle", "wb") as f: pickle.dump(data, f, protocol = 4)


#%%


f_tree = dir_halo + "tree_py.dat"

dir_fig = base + "figs/"
# Open the file and call pickle.dump.
#with open("f.pickle", "wb") as f: pickle.dump(data, f, protocol = 4)
# protocol 4 is added in python3.400
# protocol = 2 -> binary writing taking up less space. > python2.3

#data = None
# Open the file and call pickle.load.
f_halo = dir_halo + str(nout)
with open(f_halo, "rb") as f_halo:
	data = pickle.load(f_halo)
	
gals = tru.gal_list(data)			
#%%

# Make galaxy 

				