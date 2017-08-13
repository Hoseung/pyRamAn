import numpy as np
import pickle
from glob import glob

def idxfromid(org_ids, ids_ref, idxs_ref):
    i_sort = np.argsort(ids_ref)
    ids_ref_new = ids_ref[i_sort]
    idxs_ref_new = idxs_ref[i_sort]

    return idxs_ref_new[np.argsort(org_ids).argsort()]


def get_all_results(nouts,
                    prg_dir = "./test_direct_prgs_gal/",
                    out_dir = "./lambda_results/",
                    fix_idx_nout=9999):
    # ALL ALL results.
    all_sample_ids=pickle.load(open(prg_dir + "all_sample_ids.pickle", "rb"))
    all_sample_idxs=pickle.load(open(prg_dir + "all_sample_idxs.pickle", "rb"))
    allresults=[]
    for nout in nouts:
        allresults_thisnout = []
        # Better that the list of sample is stored in a separate file, instead of
        # trying to read ALL files in a directory...
        fn_all = glob(out_dir+"{}/result_sub_sample_{}_*.pickle".format(nout, nout))
        # Get right IDx
        for fn in fn_all:
            idsnow = np.array(all_sample_ids[str(nout)])
            idxsnow = np.array(all_sample_idxs[str(nout)])

            # Some results have right idx, some are wrong...
            this_result = pickle.load(open(fn, "rb"))
            allidxs = np.array([agal.idx for agal in this_result])
            if max(allidxs) < 1e6:
                allidxs = idxsnow[mtc.match_list_ind(idsnow, allidxs)]
                #idxfromid(allidxs, idsnow, idxsnow)
                for idx, agal in zip(allidxs, this_result):
                    agal.idx =idx
            allresults_thisnout.extend(this_result)
            # Don't like using glob.
            # file name will be "ixiyiz" instead of "from xxxx".
            # then just i=0:999
            #if nout >= fix_idx_nout:
            # Fix IDx on the fly.

            #allidxs = idxfromid(np.array([agal.idx for agal in allresults_thisnout]), idsnow, idxsnow)
            #for idx, agal in zip(allidxs, allresults_thisnout):
            #    agal.idx =idx
        allresults.append(allresults_thisnout)
    return allresults
