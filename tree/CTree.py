# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 17:10:10 2015

@author: hoseung
"""
from tree import treemodule
from tree import treeutils
import numpy as np
#%%
trees = treemodule.CTree()
wdir = '/home/hoseung/Work/data/AGN2/'
trees.load(filename= wdir + 'rhalo/tree_0_0_0.dat')
#%%
prgs, inds = treeutils.get_main_prg(trees, haloind=11)
import matplotlib.pyplot as plt
#%%
import tree.treeplots as trp
a = trees.data[inds]
trp.plot_all(a, a['id'][0], save=True, out_dir=wdir)

#%%
class tree_summary():
    import numpy as np
    def __init__(self, full_tree):
        self.roots = self.get_root_halo_id(full_tree)
        data = np.recarray(len(self.roots,), 
                dtype=[('ID', int), ('aexp_ini', float), ('aexp_fin', float), ('N_merger', int)])
        data['ID'] = self.roots
        #data['aexp_ini']
    
    def get_root_halo_id(self, full_tree):
        return np.unique(full_tree['tree_root_id'])
    