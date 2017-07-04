import pickle
import numpy as np
from utils import hagn
from general import defaults
import tree.halomodule as hmo


def load_idx_list(wdir="./", test=False):
    if test:
        return np.genfromtxt(wdir+"test_direct_prgs_gal/final_idxs_allmassive_gal.txt",
                 dtype=[('idx', 'i8'), ('id', 'i8')])
    else:
        return np.genfromtxt(wdir+"all_direct_prgs_gal/final_idxs_allmassive_gal.txt", dtype=int)


def all_id_to_dict(all_idxs, nnza, dir_ids = "pidxs/"):
    all_pidxs = [[] for i in range(nnza.nnza["nstep"].max())]
    all_pids = [[] for i in range(nnza.nnza["nstep"].max())]
    for final_idx in all_idxs:
        pidxs_this = pickle.load(open(dir_ids+"IDx/IDxall_direct_prgs_" + str(final_idx) +"_gal.pickle", "rb"))
        for i, pt in enumerate(pidxs_this):
            if len(pt) > 0:
                all_pidxs[i].extend(pt)
        pids_this = pickle.load(open(dir_ids+"./ID/IDall_direct_prgs_" + str(final_idx) +"_gal.pickle", "rb"))
        for i, pt in enumerate(pids_this):
            if len(pt) > 0:
                all_pids[i].extend(pt)

    all_ids=dict()
    for ap, nouts in zip(all_pids, nnza.nnza["nout"]):
        all_ids.update({str(nouts):ap})
    all_idxs=dict()
    for ap, nouts in zip(all_pidxs, nnza.nnza["nout"]):
        all_idxs.update({str(nouts):ap})
    
    return all_ids, all_idxs

