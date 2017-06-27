#import matplotlib
#matplotlib.use("Qt5Agg")
#import matplotlib.pyplot as plt
import load 
import numpy as np
import tree.halomodule as hmo
from general import defaults
from galaxymodule import rotation_parameter, rd_GM
from utils import hagn
from multiprocessing import Pool
from rot2 import cell_chunk_module as ccm
import pickle
#from rot2 import prg_modules as prm

if __name__ == "__main__":

    # for nout in nouts:
    nout = 782 
    gcat = hmo.Halo(nout=nout, is_gal=True)

    prg_dir = "./all_direct_prgs_gal/"
    all_sample_ids=pickle.load(open(prg_dir + "all_sample_ids.pickle", "rb"))
    ccm.cat_only_relevant_gals(gcat, all_sample_ids, nout) 

    args =[]
    for cat_chunk in ccm.domain_decompose_cat(gcat, nbins=5):
        args.append([cat_chunk, nout])

    # Parallelize
    with Pool(processes=4) as pool:
        pool.starmap_async(ccm.do_work, args)

    print("DONE")

    #for sub_sample in domain_decompose_cat(gcat, nbins=5)[:1]:
    #    do_work(sub_sample, nout)


