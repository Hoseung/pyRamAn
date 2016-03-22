# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 22:43:08 2015

@author: hoseung
"""

import numpy as np
import tree.treeutils as tru
import tree.treeplots as trp
import tree.rshalo as rsh
import utils.sampling as smp
from utils import util
import pickle


# Load tree
wdir = '/home/hoseung/Work/data/AGN2/'
dir_halo = wdir + "rhalo/rockstar_halos/"
f_tree = wdir + "rhalo/tree.pickle"
fn_halo = dir_halo + 'halos_81.ascii'

nout_ini = 0
nout_fi = 81
nouts = range(nout_fi, nout_ini, -1)
Nnouts = len(nouts)

# RS tree
with open(f_tree, "rb") as ft:
    rstree = pickle.load(ft)



#%%