#from rot2.prg_modules import *
import numpy as np
import pickle
import tree.halomodule as hmo
from glob import glob
from utils import hagn

def all_id_to_dict(nnza, fid_list=None, adp_dir="./all_trees/"):
    nsteps=[]
    ids=[]
    idxs=[]
    #if fid_list is None:
    all_files = glob(adp_dir+"*_adp.pickle")
    for fn in all_files:
        adp = pickle.load(open(fn, "rb"))
        if fid_list is not None:
            print(adp[0][5]["id"])
            if not adp[0][5]["id"] in fid_list:
                continue
        for this_sats in adp:
            for sat in this_sats:
                nsteps.extend(sat["nstep"])
                ids.extend(sat["id"])
                idxs.extend(sat["idx"])

    print(ids)
    print(idxs)
    print(len(ids), len(idxs))
    data = np.column_stack((nsteps, ids, idxs))#, dtype=[("nstep", "<i8"), ("id", "<i8"),("idx", "<i8")])
    data_sorted = data[np.argsort(data[:,0])]

    count, bins = np.histogram(data_sorted[:,0], bins=np.arange(data_sorted[-1,0]+2)) # 
    cum_cnt = np.cumsum(count)

    all_ids=dict()
    all_idxs=dict()
    # len(cum_cnt) + 2 = len(bins)
    #print(bins[1:-1])
    #print(cum_cnt[1:]-cum_cnt[:-1])
    for icc, nstep in enumerate(bins[1:-1]):
        nout = nnza.step2out(nstep)
        #print(nout)
        all_ids.update({str(nout):data_sorted[cum_cnt[icc]:cum_cnt[icc+1],1]})
        all_idxs.update({str(nout):data_sorted[cum_cnt[icc]:cum_cnt[icc+1],2]})

    return all_ids, all_idxs


if __name__=="__main__":
    """
        Convert per tree prgs into per snapshot galaxy id list.
    """
    nnza = hagn.Nnza()#fname="./GalaxyMaker/gal/nout_nstep_zred_aexp.txt")
    # All final IDX
    #all_idxs = load_idx_list(test=test)["idx"]

    prg_dir = "./all_fine_direct_prgs_gal/"
    #    prg_dir = "./all_direct_prgs_gal/"
    # Load all IDx/IDxall_direct_prgs (and IDs)
    # and gather per each step.
    #all_sample_ids, all_sample_idxs=all_id_to_dict(all_idxs, nnza, prg_dir + "pidxs/")
    etg_fid = np.genfromtxt("RedGals.txt", dtype=int)
    all_sample_ids, all_sample_idxs=all_id_to_dict(nnza, fid_list = etg_fid, adp_dir=prg_dir)
    pickle.dump(all_sample_ids, open(prg_dir + "etg_sample_ids.pickle","wb"))
    pickle.dump(all_sample_idxs, open(prg_dir + "etg_sample_idxs.pickle","wb"))

