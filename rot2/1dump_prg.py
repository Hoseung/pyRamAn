import pickle
import numpy as np
import tree
import os
from rot import refined_tree

def mkdir(dirpath):
    if not os.path.isdir(dirpath):
        os.mkdir(dirpath)

def fname(fidx, is_gal=False):
	if is_gal:
		return "all_direct_prgs_"+str(fidx)+"_gal.pickle"
	else:
		return "all_direct_prgs_"+str(fidx)+"_hal.pickle"

if __name__ == "__main__":

    test = True

    wdir = './'
    is_gal = True
    # Do I ever care halos?

    if is_gal:
        Mcut = 1e10
        if test:
            basedir = wdir + "test_direct_prgs_gal/"
        else:
            basedir = wdir + "all_direct_prgs_gal/"
    if not is_gal:
        Mcut = 5e10
        if test:
            basedir = wdir + "test_direct_prgs_gal/"
        else:
            basedir = wdir + "all_direct_prgs_hal/"

    mkdir(basedir)
    #mkdir(basedir+"pidxs/")
    #mkdir(basedir+"pidxs/ID/")
    #mkdir(basedir+"pidxs/IDx/")
    #mkdir(basedir+"ptrees/")

    # Load data

    tt = tree.tmtree.Tree(is_gal=is_gal)
    print("Loading tree done")

    tnow = tt.tree[tt.tree["nstep"]==max(tt.tree["nstep"])]
    large_last = tnow[tnow["m"] > Mcut]# * (tnow["m"] < 2)]  above 1e10

    # TEST
    if test:
        large_last = large_last[(np.abs(large_last["xp"][:,0])
                               + np.abs(large_last["xp"][:,1])
                               + np.abs(large_last["xp"][:,2]) < 100)] # 

    print("Number of sample galaxies:", len(large_last))

    final_idxs = large_last["idx"]
    final_ids = large_last["id"]

    np.savetxt(basedir + "final_idxs_allmassive_gal.txt", np.c_[final_idxs,final_ids], fmt='%d  %d')

    num_gal = len(final_ids)

    for i, (fid,fidx) in enumerate(zip(final_ids, final_idxs)):
        #maintree, idx_prgs_alltime, id_prgs_alltime = tt.extract_direct_full_tree(fidx, return_id=True)
        # All list
        #pickle.dump(idx_prgs_alltime, open(basedir+"pidxs/IDx/IDx"+fname(fidx, is_gal), "wb"))
        #pickle.dump(id_prgs_alltime, open(basedir+"pidxs/ID/ID"+fname(fidx, is_gal), "wb"))
        #adp = tt.get_all_trees(idx_prgs_alltime)
        #pickle.dump((maintree, idx_prgs_alltime, adp), open(basedir+"ptrees/"+fname(fidx, is_gal), "wb"))
        refined_tree.refined_tree(tt, fidx,
                                  f_dist_sum = 0.66
                                  r_fi = 0.5
                                  step_early_enough = 30
                                  m_small_enough = 3.3e8
                                  too_short_ref = 2
                                  threshold_score=2.0
                                  l_tree_max=50
                                  l_too_short_close=15,
                                  do_plot=True,
                                  out_dir=basedir)
        print("{}-th / {}".format(i, num_gal), end="\r")
