#import matplotlib
#matplotlib.use("Qt5Agg")
#import matplotlib.pyplot as plt
#import load 
import numpy as np
import tree.halomodule as hmo
#from utils import hagn
from multiprocessing import Pool
from rot2 import cell_chunk_module as ccm
import pickle
#from rot2 import prg_modules as prm

if __name__ == "__main__":

    nnza = np.genfromtxt("./nout_nstep_zred_aexp.txt",
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])

    # for nout in nouts:
    nouts = nnza["nout"]
    nout = nouts[1]
    gcat = hmo.Halo(nout=nout, is_gal=True)

    prg_dir = "./all_direct_prgs_gal/"
    all_sample_ids=pickle.load(open(prg_dir + "all_sample_ids.pickle", "rb"))
    ccm.cat_only_relevant_gals(gcat, all_sample_ids, nout) 

    args =[]
    for i, cat_chunk in enumerate(ccm.domain_decompose_cat(gcat, nbins=10)):
        args.append([cat_chunk, nout])
    #    ccm.do_work(cat_chunk, nout)

    # Parallelize
    with Pool(processes=3) as pool:
        pool.starmap(ccm.do_work, args)
        # Why can't I use starmap_async?

    print("DONE")

    #for sub_sample in domain_decompose_cat(gcat, nbins=5)[:1]:
    #    do_work(sub_sample, nout)


