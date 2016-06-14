

from analysis.cal_lambda import *
import pickle
import numpy as np
import load


# load prg_only_tree
wdir = './'
prg_only_tree = pickle.load(open(wdir + "prg_only_tree.pickle", 'rb'))

#
nouts = range(187, 180, -1)
for nout in nouts:
    info = load.info.Info(nout=nout, base=wdir, load=True)
    #if nout < 187:
    mstar_min = 2 * get_mstar_min(info.aexp)

    allgal, allhal = get_sample_gal(wdir, nout, info, prg_only_tree, mstar_min)

    print(allgal, allhal) 
