#from rot2.prg_modules import *
import numpy as np
import pickle
import tree.halomodule as hmo
from glob import glob

def all_id_to_dict(adp_dir="./all_trees/"):
    nsteps=[]
    ids=[]
    idxs=[]
    all_files = glob(adp_dir+"*_adp.pickle")
    for fn in all_files:
        maintree, idx_prgs_alltime, id_prgs_alltime, adp = pickle.load(open(fn, "rb"))
        for this_sats in adp:
            for sat in this_sats:
                nsteps.extend(sat["nstep"])
                ids.extend(sat["id"])
                idxs.extend(sat["idx"])

    data = np.column_stack((nsteps, ids, idxs))#, dtype=[("nstep", "<i8"), ("id", "<i8"),("idx", "<i8")])
    data_sorted = data[np.argsort(data[:,0])]

    count, bins = np.histogram(data_sorted[:,0], bins=np.arange(data_sorted[-1,0]+2)) # 
    cum_cnt = np.cumsum(count)

    all_ids=dict()
    all_idxs=dict()
    # len(cum_cnt) + 2 = len(bins)
    for icc, nout in enumerate(bins[1:-1]):
        all_ids.update({str(nout):data_sorted[cum_cnt[icc]:cum_cnt[icc+1],1]})
        all_idxs.update({str(nout):data_sorted[cum_cnt[icc]:cum_cnt[icc+1],2]})

    return all_ids, all_idxs


if __name__=="__main__":
    """
        Convert per tree prgs into per snapshot galaxy id list.
    """
    test = True
    nnza = hagn.Nnza()
    # All final IDX
    #all_idxs = load_idx_list(test=test)["idx"]

    if test:
        prg_dir = "./test_direct_prgs_gal/"
    else:
        prg_dir = "./all_direct_prgs_gal/"
    # Load all IDx/IDxall_direct_prgs (and IDs)
    # and gather per each step.
    #all_sample_ids, all_sample_idxs=all_id_to_dict(all_idxs, nnza, prg_dir + "pidxs/")
    all_sample_ids, all_sample_idxs=all_id_to_dict(apd_dir="./all_trees/")
    pickle.dump(all_sample_ids, open(prg_dir + "all_sample_ids.pickle","wb"))
    pickle.dump(all_sample_idxs, open(prg_dir + "all_sample_idxs.pickle","wb"))

