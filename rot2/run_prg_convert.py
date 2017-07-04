from rot2.prg_modules import *
from general import defaults
import numpy as np
import pickle
import tree.halomodule as hmo


if __name__=="__main__":
    """
        Convert per tree prgs into per snapshot galaxy id list.
    """
    test = True
    dfl = defaults.Default()
    nnza = hagn.Nnza()

    all_idxs = load_idx_list(test=test)["idx"]

    if test:
        prg_dir = "./test_direct_prgs_gal/"
    else:
        prg_dir = "./all_direct_prgs_gal/"

    all_sample_ids, all_sample_idxs=all_id_to_dict(all_idxs, nnza, prg_dir + "pidxs/")
    pickle.dump(all_sample_ids, open(prg_dir + "all_sample_ids.pickle","wb"))
    pickle.dump(all_sample_idxs, open(prg_dir + "all_sample_idxs.pickle","wb"))

    nout = 782
    gcat = hmo.Halo(nout=782, is_gal=True)
    subsample_catalog(gcat, all_sample_ids, nout)

