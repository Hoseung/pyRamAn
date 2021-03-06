import numpy as np
import tree.halomodule as hmo
from multiprocessing import Pool
from rot2 import cell_chunk_module as ccm
import pickle
import os
import utils.match as mtc

if __name__ == "__main__":
    test=False
    nnza = np.genfromtxt("./nout_nstep_zred_aexp.txt",
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])

    nouts = nnza["nout"]
    
    # nouts with Hydro raw data.
    nnza_cell = np.genfromtxt("./nout_nstep_zred_aexp_63.txt",
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])

    nbins=10
    save_cell=True
    wdir = '/home/hoseung/data/Horizon-AGN/'
    out_base='./'
    if test:
        prg_dir = "./test_fine_direct_prgs_gal/"
    else:
        prg_dir = "./all_fine_direct_prgs_gal/"
 
    outdir = "./lambda_results/"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
 
    sample_ids=pickle.load(open(prg_dir + "all_sample_ids.pickle", "rb"))
    #sample_idxs=pickle.load(open(prg_dir + "all_sample_idxs.pickle", "rb"))
    #sample_ids = np.intersect1d(small_sample_id, sample_ids[str(nout)])

    #for nout in nouts:
    for nout in nnza_cell["nout"][:13]:
        gcat = hmo.Halo(nout=nout, is_gal=True)
        if not os.path.isdir(outdir + str(nout)):
            os.mkdir(outdir + str(nout))
        if save_cell:
            save_cell_dir=out_base+"CELL_"+str(nout).zfill(5)
            if not os.path.isdir(save_cell_dir):
                print("mkdir", save_cell_dir)
                os.mkdir(save_cell_dir)
         
        #ccm.cat_only_relevant_gals(gcat, sample_ids, sample_idxs, nout)
        allgal_now = np.unique(sample_ids[str(nout)])
        #allgal_now_idxs = np.unique(sample_idxs[str(nout)])
        #print(gcat.data.dtype)
        gcat.data = gcat.data[mtc.match_list_ind(gcat.data["id"], allgal_now)]
        #gcat.data["idx"] = allgal_now_idxs[mtc.match_list_ind(gcat.data["id"], allgal_now)]

        # Smaller test sample
        args =[]
        for i, cat_chunk in enumerate(ccm.domain_decompose_cat(gcat, nbins=nbins)):
            if i < 0:
                continue
            if len(cat_chunk) == 0:
                print("No target galaxy in this sub volume")
                continue
            args.append([cat_chunk, nout, i])
            #print(i,cat_chunk[0]["id"])
            #if test:
            ccm.do_work(cat_chunk, nout, i, do_cell=False, base=wdir)
 
        print("# of subsamples", len(args))
        # Parallelize
        if True:
            with Pool(processes=4) as pool:
                pool.starmap(ccm.do_work, args)
            # Why can't I use starmap_async?
 
        print("DONE")
 
        #for sub_sample in domain_decompose_cat(gcat, nbins=5)[:1]:
        #    do_work(sub_sample, nout)
