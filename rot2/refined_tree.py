import matplotlib.pyplot as plt
#import tree
import pickle
import numpy as np
from rot2 import treetmp
from rot2 import fix_tree

#H0=70.4
#tt= tree.tmtree.Tree(is_gal=True)
#treetmp.build_kdtree(tt, H0)
#tt = pickle.load(open("tree_kd.pickle","rb"))

#tnow = tt.tree[tt.tree["nstep"] == tt.tree["nstep"].max()]
#large_gals = tnow["idx"][tnow["m"] > 1e10]
# FYI, gal ID is compatible with the coarse tree!

# Fix tree parameters

def refined_tree(tt, fidx, f_dist_sum = 0.66,
                          r_fi = 0.5,
                          step_early_enough = 30,
                          m_small_enough = 3.3e8,
                          too_short_ref = 2,
                          threshold_score=2.0,
                          l_tree_max=50,
                          l_too_short_close=15,
                          do_plot=False,
                          out_dir="./"):
    """
    Parameters
    ----------
        f_dist_sum :
            a satellite inside the f_dist_sum*(r1+r2) is considered to have merged.
        r_fi :
            dd
        step_early_enough : 
            A tree starting earlier than this is thought to be complete.
        m_small_enough : 
            smaller than this mass is too small that I can ignore.
        too_short_ref : 
            A tree let shorter than this is not viable for fitting. Ignored.
        threshold_score : 
            Threshold for test_similar.
        l_tree_max : 
            pre-merger state of the satellite is assumed to be < l_tree_max steps back. 
        l_too_short_close :
            If a tree is shorter than this and always within f_dist_sum, 
            this is a remnant of previous merger (which is disconnected from 
            this sat tree, and likely exists separately.) Because the main body of the sat
            ,including pre-merger state, presumably presents, this remnant is of no use. remove.
        do_plot :
            plot trees before and after fixing. 
        
    """

    max_dM_frac1=100.0
    m_frac_min1=0.3
    max_dM_frac2=10.0
    m_frac_min2=0.3

    fix_result = []
    print("This FIDX", fidx)
    
    # gather idxs of all sat roots(final IDX).
    #   # 100 : [[ [main] ],
    #   #  99 :  [ [idx_sat1], [idx_sat2], ... ],
    #   #  98 :  [ [idx_sat1], [idx_sat2], ... ],
    #            ...]
    maintree, idx_prgs_alltime, id_prgs_alltime, macc \
        = treetmp.extract_direct_full_tree(tt, fidx,
                                   return_id=True,
                                   return_macc=True,
                                   max_dM_frac=max_dM_frac1,
                                   m_frac_min = m_frac_min1)
    
    # Based on idx lists, 
    # get main trees of all sat roots.
    # tree = np.ndarray.
    #        dtype=["nstep","xp","vp","lp","cvel",....] # ~20 fields?
    #
    #   # 100 : [[ [tree_main] ],
    #   #  99 :  [ [tree_sat1], [tree_sat2], ... ],
    #   #  98 :  [ [tree_sat1], [tree_sat2], ... ],
    #            ...]
    #
    # Redundant sat idxs are removed. 
    # (when a idx at previous step is found in sat trees that start from earlier steps)
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #               |
    #               |
    #               |
    #
    #
    #
    #
    #
    #
    #  Parameters:
    #      max_dM_frac : Less than by this factor a halo can change in mass. 
    #                    20-fold mass increase/decrease is unphysical, of course.
    #                    But underestimatig halo mass during merger can result in 
    #                    that much change. (Usually up to ~5x in galaxy merger trees)
    #      m_frac_min : More than this fraction of a halo in the previous step must be
    #                   transferred to the descendant to be the main progenitor. 
    #                   If it's binary merger, m_frac_min must be >0.5
    #                   But there could be multiple mergers. 
    adp = treetmp.get_all_trees(tt,
                                idx_prgs_alltime,
                                verbose=False,
                                max_dM_frac=max_dM_frac2,
                                m_frac_min =m_frac_min2)
    adp[0]=[maintree]
    if do_plot:
        treetmp.check_tree(adp, save=True, nstep_min = 50, figure_type="simple")
    #treetmp.check_tree(adp, save=True, figure_type="detailed")
 
    # FIX
    for k, sats_now in enumerate(adp):
        if k ==0:
        #if k < 322:
            continue
        len_sats = len(sats_now)
        if len_sats > 0:
            #print("This step", sats_now[0][0]["nstep"])
            j=0
            j_skip=0
            while j + j_skip < len_sats:
                sat = sats_now[j]
                if len(sat) <= too_short_ref:
                    adp[k].pop(j)
                    # One element is taken from the list.
                    j_skip +=1
                    continue
                elif sat[-1]["m"] < m_small_enough and sat[-1]["nstep"] < 2*step_early_enough:
                    # Keep
                    j+=1
                    continue
                elif sat[-1]["nstep"] < step_early_enough:
                    # Keep
                    j+=1
                    continue
                else:
                    all_dist = fix_tree.main_sat_dist(maintree, sat, ratio=True)
                    merger_fi = np.argmax(all_dist > r_fi)
                    #print("\n merger_fi", merger_fi)
 
                    l_org = len(sat)
                    if all_dist[0] < f_dist_sum and merger_fi==0:
                        print("Ignore", all_dist[:10])
                        # Always closer than the thresold.
                        # Properties of the satellite should have been measured earlier.
                        # Only the final step of the tree is useful to track
                        # the end of the merger.
                        sat["nextsub"] = -123
 
                    # Maybe not long enough.
                    # Although 50 is very conservative value. (~2Gyr)
                    # broken tree!
                    print("\n Fixing a broken tree", sat["nstep"][0], sat["idx"][0], k,j)
                    fix_result.append(fix_tree.fix_broken_tree(adp,k,j,
                                      fields=["m", "rho_0", "cvel", "spin", "ek", "rs"],
                                      dist_tol=0.5,
                                      dnstep_max_connect=7,
                                      nstep_max=maintree[0]["nstep"],
                                      threshold_score=threshold_score))
 
                    #adp[k][j] = adp[k][j][:l_tree_max]
                    if fix_result[-1] != "good":
                        j +=1
                    # If fix_broken_tree was successful,
                    # the tree has been modified and the end of the tree
                    # is now to be merged with another, earlier treelet.
                    # So, do not increment j and let it iterate again.
 
    #treetmp.check_tree(adp, save=True,
    #                   nstep_min = 50,
    #                    figure_type="simple",
    #                    suffix="middle")

    # remove small treelets 
    for k, sats_now in enumerate(adp):
        if k ==0:
            continue
        len_sats = len(sats_now)
        if len_sats > 0:
            j=0
            j_skip=0
            while j + j_skip < len_sats:
                sat = sats_now[j]
                if len(sat) <= too_short_ref:
                    adp[k].pop(j)
                    j_skip +=1
                    continue
                elif sat[-1]["m"] < m_small_enough or sat[-1]["nstep"] < step_early_enough:
                    j+=1
                    continue
                elif sat[-1]["nstep"] < step_early_enough:
                    j+=1
                    continue
                else:
                    all_dist = fix_tree.main_sat_dist(maintree, sat, ratio=True)
                    merger_fi = np.argmax(all_dist > r_fi)
                    #print("merger_fi", merger_fi)
 
                    l_org = len(sat)
                    if all_dist[0] < f_dist_sum and merger_fi==0 and len(sat) < l_too_short_close:
                        print("Ignore", all_dist[:10])
                        # Always closer than the thresold.
                        # Properties of the satellite should have been measured earlier.
                        # Only the final step of the tree is useful to track
                        # the end of the merger.
                        adp[k].remove(sat)
                        j_skip+=1
                    else:
                        j+=1
 
    if do_plot:
        treetmp.check_tree(adp,
                       save=True,
                        suffix="fixed",
                        nstep_min = 50,
                        figure_type="simple")
        plt.close()
 
    pickle.dump(adp,open(out_dir+"{}_adp.pickle".format(maintree[0]["idx"]), "wb"))


    # How many are "no candidates", "no good matches", or "good"?
    #print(fix_result.count("no candidates"),
    #    fix_result.count("no good matches"),
    #    fix_result.count("good"))
