# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 14:32:50 2015

generate main progenitor list of a given galaxy

@author: hoseung
"""


import multiprocessing as mp    
import load
from tree import tmtree
import numpy as np
import utils.sampling as smp
import tree.halomodule as hmo 
from utils.util import reimport

wdir = '/home/hoseung/Work/data/AGN2/'
ftree = wdir + 'rhalo/tree.pickle'
nout_ini = 37
nout_fi = 187

#%% load tree data
#import pickle
#f = open(ftree, "rb")
#tt = pickle.load(f)

#%% or,
import tree.treemodule as tmo
reimport(tmo)
tt = tmo.CTree()
tt.load(ftree)
# If .pickle file is given,
# tree  data is loaded not from ascii file, but from pickle.
#%%
import tree.treeutils as tru
prgs = tru.get_main_prg(tt, haloid=1634975) # id? org_id? what...??