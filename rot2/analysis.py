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


def smooth(x, beta=5, window_len=20, monotonic=False, clip_tail_zeros=True):
    """
    kaiser window smoothing.

    If len(x) < window_len, window_len is overwritten to be len(x).
    This ensures to return valid length fo array, but with modified window size.

    Parameters
    ----------
        window_len = 20

        monotoinc =

        clip_tail_zereos = True
            returned array is shorter than the original.

    beta = 5 : Similar to Hamming


    """
    lx = len(x)
    if clip_tail_zeros:
        x = x[:max(np.where(x > 0)[0])+1]

    if monotonic:
        """
        if there is an overall slope, smoothing may result in offset.
        compensate for that.
        """
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y=np.arange(len(x)))
        xx = np.arange(len(x)) * slope + intercept
        x = x - xx

    # extending the data at beginning and at the end
    # to apply the window at the borders
    window_len = min([window_len, len(x)])
    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]] # concatenate along 0-th axis.
    # periodic boundary.
    w = np.kaiser(window_len,beta)
    y = np.convolve(w/w.sum(), s, mode='valid')
    aa = len(y)-lx
    if monotonic:
        return y[int(window_len/2):len(y)-int(window_len/2) + 1] + xx
    else:
        #return y[int(window_len-1/2)-2:len(y)-int((window_len-1)/2)]
        return y[int(aa/2):int(aa/2)+len(x)]
