from prg_modules import *

if __name__=="main":
    """
        Convert per tree prgs into per snapshot galaxy id list.
    """
    dfl = defaults.Default()
    nnza = hagn.Nnza()

    all_idxs = load_idx_list()

    prg_dir = "./all_direct_prgs_gal/"

    all_sample_ids, all_sample_idxs=all_id_to_dict(all_idxs, nnza, prg_dir + "pidxs/")
    pickle.dump(all_sample_ids, open(prg_dir + "all_sample_ids.pickle","wb"))
    pickle.dump(all_sample_idxs, open(prg_dir + "all_sample_idxs.pickle","wb"))

    nout = 782
    gcat = hmo.Halo(nout=782, is_gal=True)
    subsample_catalog(gcat, all_sample_ids, nout)

