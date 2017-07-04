import numpy as np
import tree.halomodule as hmo
from multiprocessing import Pool
from rot2 import cell_chunk_module as ccm
import pickle

if __name__ == "__main__":
    test=True
    nnza = np.genfromtxt("./nout_nstep_zred_aexp.txt",
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])

    # for nout in nouts:
    nouts = nnza["nout"]
    nout = nouts[0]
    gcat = hmo.Halo(nout=nout, is_gal=True)
    if test:
        prg_dir = "./test_direct_prgs_gal/"
    else:
        prg_dir = "./all_direct_prgs_gal/"

    all_sample_ids=pickle.load(open(prg_dir + "all_sample_ids.pickle", "rb"))
    if test:
        small_sample_id= gcat.data["id"][(gcat.data["x"] < 0.4) * (gcat.data["y"] < 0.4) * (gcat.data["z"] < 0.4)]
        all_sample_ids = np.intersect1d(small_sample_id, all_sample_ids)

    ccm.cat_only_relevant_gals(gcat, all_sample_ids, nout)
    # Smaller test sample

    args =[]
    for i, cat_chunk in enumerate(ccm.domain_decompose_cat(gcat, nbins=10)):
        #if i < 64:
        #    continue
        if len(cat_chunk) == 0:
            print("No target galaxy in this sub volume")
            continue
        args.append([cat_chunk, nout])
        if test:
            ccm.do_work(cat_chunk, nout)

    # Parallelize
    if not test:
        with Pool(processes=3) as pool:
            pool.starmap(ccm.do_work, args)
        # Why can't I use starmap_async?

    print("DONE")

    #for sub_sample in domain_decompose_cat(gcat, nbins=5)[:1]:
    #    do_work(sub_sample, nout)
