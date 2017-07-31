import numpy as np



def get_diff(sat1,sat2,field):
    return np.abs((sat1[field]-sat2[field])/sat1[field])+ np.abs((sat1[field]-sat2[field])/sat1[field])

def test_similar(sat1, sat2, fields):
    """
    get_diff can tell which of cancidates is the most appropriate,
    but provides no information in absolute scale.
    """
    # mass
    metrics= [get_diff(sat1,sat2,ff) for ff in fields]
    return np.sum(metrics) / len(metrics)


def interpol_treelets(ref, missing_tree,
                      ref_post,
                      fields_interpol,
                      dnstep=7, poly_deg=5):

    incomplete_data = np.concatenate((this_branch, candidate_branches[0]))
    nstep_fit_first = ref[-1]["nstep"] - dnstep
    nstep_fit_last = ref_post[0]["nstep"] + dnstep
    i_start = np.where(incomplete_data["nstep"] == nstep_fit_first)[0][0]
    i_end = np.where(incomplete_data["nstep"] == nstep_fit_last)[0][0]

    X=incomplete_data["nstep"][i_end:i_start]
    print(nstep_fit_first, nstep_fit_last)
    print(i_start, i_end)
    print(X)
    x_pred = missing_tree["nstep"]

    for this_field in fields_interpol:
        Y=incomplete_data[this_field][i_end:i_start]
        if Y.ndim==2:
            for vind in range(3):
                z = np.polyfit(X,Y[:,vind], deg=poly_deg) # Degree?
                p = np.poly1d(z)
                missing_tree[this_field][:,vind]=p(x_pred)
        else:
            z = np.polyfit(X,Y, deg=poly_deg) # Degree?
            p = np.poly1d(z)

            missing_tree[:][this_field]=p(x_pred)

def fix_broken_tree(adp,
                    istep_pre,isat_pre,
                    fields=["m", "rho_0", "cvel", "spin", "ek", "rho_c"],
                    score_good=1.0):
    """

    score = 1.0 : 50% difference in all fields.
    """

    nstep_max = 388
    dnstep_max_connect=3 # trees at maximum 5 snapshot apart can be connected.

    ref = adp[istep_pre][isat_pre]
    nstep_ini = ref["nstep"][-1]
    nstep_fi = ref["nstep"][0]
    print("Starts at {},  Ends at {}".format(nstep_ini, nstep_fi))

    # All branches ends between nstep_ini+dnstep_max_connect+1 ~ nstep_ini+1  or, 374 ~ 371
    candidate_branches=[]
    for sats_now in adp[-nstep_ini+1:-nstep_ini+1+dnstep_max_connect]:
        if len(sats_now) > 0:
            for this_sat in sats_now:
                candidate_branches.append(this_sat)
    # print(len(candidate_branches), "candidates")
    # Per step / per galaxy
    scores = []
    for cb in candidate_branches:
        scores.append(test_similiar(this_branch[-1], cb[0], fields))
            # Or, may be this_branch[-2] and cb[1] show better correlation??s

    # best match
    ibest = np.argmin(scores)
    if scores[ibest] < threshold_score:
        ref_post=candidate_branches[ibest]
    else:
        # No good enough match.
        return

    # Merge
    n_steps_missing = min(ref["nstep"][ref["nstep"]>0]) - max(ref_post["nstep"])-1
    print(n_steps_missing, min(ref["nstep"][ref["nstep"]>0]), max(ref_post["nstep"]))
    missing_tree = np.zeros(n_steps_missing, dtype=ref.dtype)
    # Fill in trivial values
    for i, tm in enumerate(missing_tree):
        tm["nstep"]=nstep_ini-i-1
        tm["id"]=-ref["id"][-1]
        tm["idx"]=-ref["idx"][-1] # if idx < -1, it is a ghost.
        tm["f_ind"]=-ref_post["idx"][0] # Not index, but refer to the original father treelet
        tm["s_ind"]-ref["idx"][-1] # Not index, but refer to the original son treelet
        tm["nsons"]=1
        tm["nprgs"]=1
    # Interpolate Missing steps
    fields_interpol = ["m", "xp", "vp", "lp", "ek", "ep", "et", "rho_c", "rho_0"]
    interpol_treelets(ref, missing_tree, ref_post, fields_interpol, dnstep=7, poly_deg=5)
    New_tree = np.concatenate((ref, missing_tree, ref_post))
    # Finally, replace the original tree
    adp[istep_pre][isat_pre] = New_tree

    print(-ref_post[0]["nstep"]-1)
    adp[-ref_post[0]["nstep"]-1].remove(ref_post)
    #ref_post = None
    return ref_post


def main_sat_dist(main, sat, ratio=True):
    main_part = main[np.in1d(main["nstep"], sat["nstep"])]

    d_pos = sat["xp"]-main_part["xp"]
    if ratio:
        return np.sqrt(np.square(sat["xp"][:,0]-main_part["xp"][:,0])
                  +np.square(sat["xp"][:,1]-main_part["xp"][:,1])
                  +np.square(sat["xp"][:,2]-main_part["xp"][:,2]))/(main_part["rvir"]+sat["rvir"])
    else:
        return np.sqrt(np.square(sat["xp"][:,0]-main_part["xp"][:,0])
                  +np.square(sat["xp"][:,1]-main_part["xp"][:,1])
                  +np.square(sat["xp"][:,2]-main_part["xp"][:,2]))

def ind_pericentric_dist_ok(main, sat, threshold_ratio=2.0):
    """
    Return the index of the last perigee where d/(r1+r2) > threshold_ratio,
    Or false.
    """
    all_dist = main_sat_dist(main, sat, ratio=True)
    # find local minima
    #d_min = argrelmax(all_dist, mode="clip")[-1]
    return all_dist[-1]>threshold_ratio
# Automatically  detect broken trees and attempt to fix them.

# It is sufficient that a tree is complete up to the point
# when the "pericentric" distance of the sat and main becomes smaller than
# n times the sum of the radii of the two.

f_dist_sum = 2.5
step_early_enough = 30
m_small_enough = 3.3e8

for i, sats_now in enumerate(adp):
    if i ==0:
        continue
    if len(sats_now) > 0:
        print("This step", sats_now[0][0]["nstep"])
        for j, sat in enumerate(sats_now):
            dist_ok=ind_pericentric_dist_ok(main, sat, threshold_ratio=f_dist_sum)
            if (sat["nstep"][-1] > step_early_enough) and (sat["m"][-1] > m_small_enough):
                dist_ok=ind_pericentric_dist_ok(main, sat, threshold_ratio=f_dist_sum)
                if not(dist_ok):
                    # broken tree!
                    print("Broken", sat["nstep"][0],sat["idx"][0], sat["idx"][-1])
                    fix_broken_tree(adp,i,j,fields=["m", "rho_0", "cvel", "spin", "ek", "rho_c"])
