# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 10:57:51 2015

@author: hoseung
"""

# load halo
import tree.load_halo as h
import tree.treeutils as tru
import tree.treeplots as trp
import pickle
hbase = '/home/hoseung/Work/data/AGN2/halos_py/'
with open(hbase + 'halos_000.pickle', 'rb') as f:
    halo = pickle.load(f)


class Tree:
    def __init__(self):
        pass
    
    def import_data(self,)


# What whould you like to have in class?

# lengh of 















from pylab import *

plot(halo['id'])
show()

#%%


X