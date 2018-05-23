import numpy as np
import matplotlib.pyplot as plt
import utils.match as mtc

def get_diff_0(sat1,sat2,field):
    return np.abs((sat1[field]-sat2[field])/sat1[field])+ np.abs((sat1[field]-sat2[field])/sat2[field])

def get_diff(sat1,sat2,field):
    return np.abs((sat1[field]-sat2[field])/sat1[field])+ np.abs((sat1[field]-sat2[field])/sat2[field])


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
                      dnstep=7, poly_deg=3):

    incomplete_data = np.concatenate((ref, ref_post))
    nstep_fit_first = max([ref_post[-1]["nstep"],ref[-1]["nstep"] - dnstep])
    nstep_fit_last = min([ref[0]["nstep"],ref_post[0]["nstep"] + dnstep])
    i_start = np.where(incomplete_data["nstep"] == nstep_fit_first)[0][0]
    i_end = np.where(incomplete_data["nstep"] == nstep_fit_last)[0][0]

    X=incomplete_data["nstep"][i_end:i_start]
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

def interpol_treelets_multi(ref, missing_tree,
                      refs_post,
                      fields_interpol,
                      dnstep=7, poly_deg=3):
    """
    ref_post is a list of treelets.
    """
    isort_post_nstep_last = np.argsort([cb["nstep"][0] for cb in refs_post])
    ref_ordered = [refs_post[i] for i in isort_post_nstep_last[::-1]]
    # list of treelets. later comes first.
    ref_ordered.insert(0, ref)

    #incomplete_data = np.concatenate(ref_ordered)
    n_missing_first = missing_tree["nstep"][0]
    for i in range(len(ref_ordered)-1):
        t_l = ref_ordered[i]
        t_r = ref_ordered[i+1]

        if t_l["nstep"][-1]-1 < t_r["nstep"][0]:
            continue
        incomplete_data = np.concatenate(ref_ordered[i:i+2])

        nstep_fit_first = max([t_r[-1]["nstep"],t_l[-1]["nstep"] - dnstep])
        nstep_fit_last = min([t_l[0]["nstep"],t_r[0]["nstep"] + dnstep])
        #i_start = np.where(incomplete_data["nstep"] == nstep_fit_first)[0][0]
        #i_end = np.where(incomplete_data["nstep"] == nstep_fit_last)[0][0]

        X=incomplete_data["nstep"][max(0,len(t_l)-dnstep):len(t_l)+dnstep]
        #print(incomplete_data["nstep"])
        #print(len(t_l),dnstep)
        #print("X", X)
        x_pred=np.arange(t_l["nstep"][-1]-1,t_r["nstep"][0],-1)
        #print("X_pred length", len(x_pred))
        i_this_t_l=n_missing_first-(t_l["nstep"][-1]-1)
        #print(" larger than 0?", i_this_t_l)
        for this_field in fields_interpol:
            Y=incomplete_data[this_field][max(0,len(t_l)-dnstep):len(t_l)+dnstep]
            #Y=np.concatenate(ref_ordered)
            if Y.ndim==2:
                for vind in range(3):
                    z = np.polyfit(X,Y[:,vind], deg=poly_deg) # Degree?
                    p = np.poly1d(z)
                    missing_tree[this_field][i_this_t_l:i_this_t_l+len(x_pred),vind]=p(x_pred)
            else:
                z = np.polyfit(X,Y, deg=poly_deg) # Degree?
                p = np.poly1d(z)
                #print("len missing_tree", len(missing_tree), x_pred)
                #print(i_this_t_l, i_this_t_l+len(x_pred))
                #print(missing_tree["nstep"],n_missing_first-(t_l["nstep"][-1]-1))
                #print(missing_tree[this_field][i_this_t_l:i_this_t_l+len(x_pred)], p(x_pred))
                #ppp =
                missing_tree[this_field][i_this_t_l:i_this_t_l+len(x_pred)]=p(x_pred)[:]

    return ref_ordered[-1]

def test_dist_similar_v(ref, cbs,
                      deg_poly=2,
                      n_pos_interpol=5,
                      nstep_back_max=5,
                      dist_tol=0.2,
                      do_plot=False):
  """
  Use velocity * dt to estimate the next position.

  Todo.
  polyfit sometimes fail.
  """
  pass


def test_dist_similar(ref, cbs,
                      deg_poly=2,
                      n_pos_interpol=5,
                      nstep_back_max=5,
                      dist_tol=0.2,
                      do_plot=False):
    """
    dist_tol : tolerance in position error. Defaults to 200Kpc
    Other properties may show better similarity once the galaxy
    recover its mass. But position is not affected by mass uderestimation, and
    thus the first step of the post tree should show the best match.
    """
    z = np.polyfit(ref["nstep"][-n_pos_interpol:],
                   ref["xp"][-n_pos_interpol:,0], deg=deg_poly)
    pxx = np.poly1d(z)
    z = np.polyfit(ref["nstep"][-n_pos_interpol:],
                   ref["xp"][-n_pos_interpol:,1], deg=deg_poly)
    pyy = np.poly1d(z)
    z = np.polyfit(ref["nstep"][-n_pos_interpol:],
                   ref["xp"][-n_pos_interpol:,2], deg=deg_poly)
    pzz = np.poly1d(z)

    X = np.arange(ref["nstep"][-1],ref["nstep"][-1]-nstep_back_max-1,-1)
    posx = pxx(X)
    posy = pyy(X)
    posz = pzz(X)

    if do_plot:
        fig, axs=plt.subplots(2,2)
        axs[0][0].plot(posx,posy, lw=2)
        axs[0][1].plot(posy,posz, lw=2)
        axs[1][0].plot(posz,posx, lw=2)
        axs[0][0].plot(ref["xp"][:,0],ref["xp"][:,1])
        axs[0][1].plot(ref["xp"][:,1],ref["xp"][:,2])
        axs[1][0].plot(ref["xp"][:,2],ref["xp"][:,0], label="ref")

    dists=[]
    #for cb in cbs:
    icb=0
    while icb < len(cbs) and len(cbs)>0:
        cb=cbs[icb]
        dist = np.sqrt(np.square(pxx(cb["nstep"])-cb["xp"][:,0])+
                       np.square(pyy(cb["nstep"])-cb["xp"][:,1])+
                       np.square(pzz(cb["nstep"])-cb["xp"][:,2]))
        if do_plot:
            axs[0][0].plot(cb["xp"][:,0],cb["xp"][:,1])
            axs[0][1].plot(cb["xp"][:,1],cb["xp"][:,2])
            axs[1][0].plot(cb["xp"][:,2],cb["xp"][:,0], label="{}".format(cb["nstep"][0]))
            axs[1][1].plot(dist)

        #print(cb[0]["nstep"], ref[-1]["nstep"], dist)
        if dist[0] > dist_tol:
            #print(cb)
            #print(cbs[3])
            cbs.pop(icb)
        else:
            #print("good cb")
            dists.append(dist[0])
            icb+=1

    if do_plot:
        axs[1][0].legend()
        plt.savefig("{}_at{}_cbs.png".format(ref[0]["idx"],ref[0]["nstep"]),dpi=300)
    return dists

def fix_broken_tree(adp,
                    istep_pre, isat_pre,
                    nstep_max = 388,
                    dnstep_max_connect=3,
                    dist_tol=0.2,
                    n_pos_interpol=5,
                    deg_poly=2,
                    poly_fit_plot=False,
                    fields=["m", "rho_0", "cvel", "spin", "ek", "rs"],
                    threshold_score=1.2):
    """

    score = 1.0 : 50% difference in all fields.

    lmax_tree : A tree does not have to be longer than 50 steps.
    """
    fields_interpol = ["m", "xp", "vp", "lp", "ek", "ep", "et", "rs", "rho_0"]
    # trees at maximum 5 snapshot apart can be connected.
    ref = adp[istep_pre][isat_pre]
    nstep_ini = ref["nstep"][-1]
    nstep_fi = ref["nstep"][0]
    print("[fix_broken_tree] Starts at {},  Ends at {}".format(nstep_ini, nstep_fi))

    # All branches ends between nstep_ini+dnstep_max_connect+1 ~ nstep_ini+1  or, 374 ~ 371
    candidate_branches=[]
    for sats_now in adp[nstep_max-nstep_ini+1:nstep_max-nstep_ini+dnstep_max_connect+1]:
        if len(sats_now) > 0:
            for this_sat in sats_now:
                #print("candidate?", this_sat["nstep"][0], this_sat["nstep"][-1], this_sat["idx"][0])
                if len(this_sat) > 1:
                    candidate_branches.append(this_sat)
                    #print(this_sat["nstep"][:3])
    if len(candidate_branches) == 0:
        return "no candidates"
    # print(len(candidate_branches), "candidates")
    # Per step / per galaxy
    scores = []
    # Measure position difference and remove galaxies too far away.
    # broken trees are usually due to heavy gravitational disturbance,
    # which entails sudden change in the motion.
    # Use only last a few points to better guess the sudden change in motion.
    # print("There are {} CBs".format(len(candidate_branches)))
    dists = test_dist_similar(ref, candidate_branches,
                              nstep_back_max=dnstep_max_connect,
                              dist_tol=dist_tol,
                              n_pos_interpol=n_pos_interpol,
                              deg_poly=deg_poly,
                              do_plot=poly_fit_plot)
    #print("Now {} CBs left".format(len(candidate_branches)))
    if len(candidate_branches) == 0:
        return "no candidates after dist cut"

    for cb in candidate_branches:
        try:
            score = test_similar(ref[-4], cb[2], fields)
        except:
            print("Fialed comparing earlier steps")
            score=1e4
        scores.append(min([test_similar(ref[-1], cb[0], fields),score]))
            # Or, may be this_branch[-2] and cb[1] show better correlation??s

    # best match
    # In desirable cases, dist match is good.
    # But the fit can go terrible during violent interactions.
    # Then do not conisder the dist match too much.

    # The following scenario is possible.
    # There are two cbs. one is long and far, the other is short and close.
    # The farther one matches best, but the closer one does good, too.
    # The shorter on fits in the gap between the ref and the farther one.
    # So, ref - shorter - longer form a tree.

    score_composite=np.array(scores)*np.array(dists)
    #print("Scores", score_composite)
    if sum(score_composite < threshold_score) > 1:
        i_comp = np.argsort(score_composite)
        ref_post = []
        nstep_post = []
        ref_post.append(candidate_branches[i_comp[0]])
        nstep_post.extend(candidate_branches[i_comp[0]]["nstep"])
        #nstep_early = best_cb[0]["nstep"]
        for icb in i_comp[1:]:
            if score_composite[icb] < threshold_score:
                if len(np.intersect1d(candidate_branches[icb]["nstep"], nstep_post))==0:
                    ref_post.append(candidate_branches[icb])
                    nstep_post.extend(candidate_branches[icb]["nstep"])

        nstep_post_min = min([cb["nstep"][0] for cb in ref_post])
        n_steps_missing = ref["nstep"][-1] - nstep_post_min-1
        if n_steps_missing > 0:
            missing_tree = np.zeros(n_steps_missing, dtype=ref.dtype)
            missing_tree["nstep"]=np.arange(ref["nstep"][-1]-1,nstep_post_min,-1)
            # Fill in the available data
            for cb in ref_post:
                ind_fit = mtc.match_list_ind(missing_tree["nstep"],cb["nstep"])
                if ind_fit is not None:
                    missing_tree[mtc.match_list_ind(missing_tree["nstep"],cb["nstep"])]=cb[:]
            for i, tm in enumerate(missing_tree):
                if tm["id"]==0:
                    tm["id"]=-ref["id"][-1]
                    tm["idx"]=-ref["idx"][-1] # if idx < -1, it is a ghost.
                    tm["f_ind"]=-ref["idx"][0] # Not index, but refer to the original father treelet
                    tm["s_ind"]-ref["idx"][-1] # Not index, but refer to the original son treelet
                    tm["nsons"]=1
                    tm["nprgs"]=1

            ref_last=interpol_treelets_multi(ref, missing_tree, ref_post, fields_interpol, dnstep=7, poly_deg=5)
            #print("ref last", ref_last["nstep"])
            New_tree = np.concatenate((ref, missing_tree, ref_last))
        else:
            # There's no room for second runner up.
            # No need to interpolate.
            New_tree = np.concatenate((ref, candidate_branches[i_comp[0]]))
        # Finally, replace the original tree
        # A tree does not have to be long.
        adp[istep_pre][isat_pre] = New_tree
        #print("REMOVE ARRAY1")
        for cb in ref_post:
            #print("Removiing...", cb[0]["nstep"], cb[-1]["nstep"])
            removearray(adp[nstep_max-cb[0]["nstep"]], cb)
    elif sum(score_composite < threshold_score) == 1:
        #print(np.argmin(score_composite))
        ref_post=candidate_branches[np.argmin(score_composite)]
        #print("Got one match")
        n_steps_missing = ref["nstep"][-1] - ref_post["nstep"][0]-1
        missing_tree = np.zeros(n_steps_missing, dtype=ref.dtype)
        for i, tm in enumerate(missing_tree):
            tm["nstep"]=nstep_ini-i-1
            tm["id"]=-ref["id"][-1]
            tm["idx"]=-ref["idx"][-1] # if idx < -1, it is a ghost.
            tm["f_ind"]=-ref["idx"][0] # Not index, but refer to the original father treelet
            tm["s_ind"]-ref["idx"][-1] # Not index, but refer to the original son treelet
            tm["nsons"]=1
            tm["nprgs"]=1
        # Interpolate Missing steps
        interpol_treelets(ref, missing_tree, ref_post, fields_interpol, dnstep=7, poly_deg=5)
        #print(" Merging ", ref["nstep"][-1], ref_post["nstep"][0])
        New_tree = np.concatenate((ref, missing_tree, ref_post))
        #print(New_tree["nstep"])
        # Finally, replace the original tree
        # A tree does not have to be long.
        adp[istep_pre][isat_pre] = New_tree
        #print("REMOVE ARRAY2")
        removearray(adp[nstep_max-ref_post[0]["nstep"]], ref_post)
    else:
        # No good enough match.
        print("No good match", scores)
        return "no good matches"

    print("GOOD\n")
    return "good"


def removearray(L,arr):
    #print("nstep cb, adp",arr["nstep"][0],L[0][0]["nstep"])
    ind = 0
    size = len(L)
    while ind != size and not L[ind][0]["idx"] == arr[0]["idx"]:
        ind += 1
    if ind != size:
        L.pop(ind)
    else:
        raise ValueError('array not found in list.')

def main_sat_dist(main, sat, ratio=True):
    main_part = main[np.in1d(main["nstep"], sat["nstep"])]
    #print(main["nstep"], main_part["nstep"], sat["nstep"])
    d_pos = sat["xp"]-main_part["xp"]
    if ratio:
        return np.sqrt(np.square(sat["xp"][:,0]-main_part["xp"][:,0])
                  +np.square(sat["xp"][:,1]-main_part["xp"][:,1])
                  +np.square(sat["xp"][:,2]-main_part["xp"][:,2]))/(main_part["rvir"]+sat["rvir"])
    else:
        return np.sqrt(np.square(sat["xp"][:,0]-main_part["xp"][:,0])
                  +np.square(sat["xp"][:,1]-main_part["xp"][:,1])
                  +np.square(sat["xp"][:,2]-main_part["xp"][:,2]))

def ind_pericentric_dist_ok(main, sat, r_ini=1.0, r_fi=1.0):
    """
    Return the index of the last perigee where d/(r1+r2) > threshold_ratio,
    Or false.
    """
    all_dist = main_sat_dist(main, sat, ratio=True)
    # First ind where d < r1+r2 = latest step
    merger_ini = np.argmax(all_dist > r_ini)
    # last ind where d > r1+r2 = earliest step

    merger_fi = np.argmin(all_dist > r_fi)
    # find local minima
    #d_min = argrelmax(all_dist, mode="clip")[-1]
    #print(all_dist)
    return merger_ini, merger_fi
