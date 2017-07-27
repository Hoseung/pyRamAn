import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import cm
from copy import copy

def build_kdtree(tt, H0):
    from scipy.spatial import cKDTree
    tt.kdts = []
    tt.ngals =[]
    for i in range(tt.nsteps):
        td = tt.tree[tt.tree["nstep"]==i+1]
        tt.kdts.append(cKDTree(td["xp"]/tt.aexps[i]*(0.01*H0)))
        tt.ngals.append(len(td))
    tt.ngals = np.array(tt.ngals)


def showtree(gal):
    fig, axs = plt.subplots(2)
    axs[0].plot(gal["m"])
    axs[1].plot(gal["vp"][:,0])
    axs[1].plot(gal["vp"][:,1])
    axs[1].plot(gal["vp"][:,2])
    plt.show()


def remove_false_prg(tt, fidxs, macc,
                     merger_mass_frac_min = 0.5,
                     verbose=False):
    """
    Any galaxy/halo that have gave more than none of mass to the son are listed as progenitors.
    But many of them are just passers-by who left a tiny fraction of their mass in the 'main' galaxy.
    For a progenitor to be a legitimate satellite, it should give most of its mass to the main galaxy (son at the next step)

    If a prg dose not transfer more than {merger_mass_frac_min} of its mass, remove it.
    """
    # remove false progenitors by checking mass transfer
    for i in np.arange(1,len(fidxs)):#zip(fidxs, fids, macc):
        idx_now = fidxs[i]
        m_son = tt.tree[fidxs[i-1][0]]["m"]
        if len(idx_now) > 1:
            if verbose:
                print("step {}, m_son {:.2e}".format(i,m_son))
            for idx, mm in zip(idx_now, macc[i]):
                mratio_sat = 0.01*mm*m_son/tt.tree[idx]["m"]
                # % of sat mass transfered to the son.
                if mratio_sat < merger_mass_frac_min:
                    print("Bad brarnch,",idx)
                    idx_now.remove(idx)

                #print(idx, mm, tt.tree[idx]["m"], ())
    return fidxs

def get_all_trees(self, idx_prgs_alltime,
                        skip_main=True,
                        filter_dup=True,
                        verbose=False,
                        **kwargs):
    """
        * For a given idx_prgs list of lists, find main progenitor tree of all entries.
        * A satellite can contribute to a host over multiple snapshots by
        given fractions of DM particles each time. In such case, the satellite
        appears in the host's progenitor tree several times.
        * Note that a 'multi-snapshot' satellite never be a main progenitor.
        However, I don't see a reason it can't be a secondary progenitor of
        another host halo. Let's just keep that in mind.

        Parameteres
        -----------
        skip_main : True
            skip main progenitor tree.

        Note
        ----
            1. About "skip_main" option.
            Main halo Tree is redundant. Main progenitor tree of the main halo at nstep = n
            includes all the main progenitors of the main halo at nstep = n-1.

        .. figure:: imgs/tmtree-get_all_trees.jpg
           :align:  center


        The main prg tree should be ignored.
        If the first gal of each snapshot is the main galaxy,
        then set skip_main =True.
        If idx_prgs_alltime is only for satellites,
        skip_main = False.
    """
    all_main_prgs=[]
    all_idxs_filter = []
    # loop over all nstep
    for j, sats_now in enumerate(idx_prgs_alltime):
        mainprgs=[]
        # satellites at each step
        for i,sat in enumerate(sats_now):
            if not skip_main or i!=0:
                #print("sat ind", i)
                if filter_dup and sat in all_idxs_filter:
                    mt = None
                    sats_now.remove(sat)
                    if verbose:
                        print("dup, remove", sat)
                else:
                    mt = extract_main_tree(self, sat, **kwargs)
                    if mt is None:
                        if verbose:
                            print("None!!", sat)
                        #print(satellite_roots)
                        sats_now.remove(sat)
                        #print(satellite_roots)
                    else:
                        #print(mt["nstep"])
                        mainprgs.append(mt)
                # update to already-member list.
                if mt is not None:
                    all_idxs_filter.extend(mt["idx"][1:])
        all_main_prgs.append(mainprgs)
        if filter_dup:
            if len(mainprgs) > int(~skip_main):
                for aa in mainprgs:
                    if len(aa) > 0 and hasattr(aa, "idx"):
                        #try:
                        all_idxs_filter.extend(aa["idx"][1:])
                        #except:
                            #print(mainprgs)
                            #return mainprgs
                    # last idx MUST remain in the prgs.
                for idxs in idx_prgs_alltime[j+1:]:
                    if len(idxs) > 1:
                        for idx in idxs[1:]:
                            if idx in all_idxs_filter:
                                idxs.remove(idx)

    return all_main_prgs


def extract_main_tree(self, idx,
                      mmin=3.3e8,
                      max_dM_frac=50.0,
                      m_frac_min=0.5,
                      verbose=False,
                      kdt_dist_upper = 0.1):
    """
    Extracts main progenitors from a TreeMaker tree.

    Criterion 1)
    A galaxy whose PEAK stellar mass is below mmin is considered to be unreliable.
    They can be a satellites, but not the trunk.
    Criterion 2)
    A galaxy can't grow by more than max_dm_frac times.
    Criterion 3)
    Likewise, a galaxy can't shrink by more than factor of ___ and remain.
    for example, 11, 11, 11, 11, 10.5, 10.2, 10.0 is possible. (log mass)
                 11, 11, 11, 11, 9.5, 9.4, 9.4 is NOT possible.

    example
    -------
    >>> tt = tmtree.Tree("tree.dat")
    >>> atree = tt.extract_main_tree(12345)

    TODO
    ----
    It works, but the try - except clause is error-prone.
    Explicitly check the end of progenitor tree and make the function more predictable.

    """
    # Tree reconstruction parameters
    nstep_back_max = 5 # At most 5 steps backwards.
    nstep_early = 10
    m_small = 5e8
     # in comoving Mpc.


    t = self.tree
    t_now = t[idx]

    # Criterion 1
    if t_now["m"] < mmin:
        if verbose:
            print("Unreliable, m < {:.2e}".format(mmin))
        return

    fatherIDx = self.fatherIDx
    fatherMass = self.fatherMass

    nstep = t_now["nstep"]
    if nstep <= 1:
        return
    nouts = [nstep]
    atree = np.zeros(nstep + 1, dtype=t.dtype)
    atree[0] = t_now

    istep_back = 0
    for i in range(1, nstep + 1):
        if istep_back > 0:
            istep_back -=1
            # skip as many as istep_back
            continue
        nstep_now = atree[i-1]["nstep"]
        idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
        i_ok = idx_father > 0
        if sum(i_ok) > 0:
            idx_father = idx_father[i_ok]
            macc_father = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]][i_ok]
            # In decending order of macc
            mass_father = np.array([t[fidx]["m"] for fidx in idx_father])
            m_frac_prg = atree[i-1]["m"] * (0.01*macc_father) / mass_father
            print("\n Father candidates before")
            [print("M_father_frac{:.2f}%  M_son_frac {:.2f}".format(100*mfrc, mac)) for mfrc, mac in zip(m_frac_prg, macc_father)]

            good_father = (m_frac_prg > m_frac_min)# * (idx_father>0)
            if len(good_father) == 0:
                idx=-2
            else:
                try:
                    macc_father = macc_father[good_father]
                    idx_father = idx_father[good_father]

                    if verbose:
                        print("\n Father candidates")
                        [print("{} {:.2f}%".format(idx, 100*mfrc)) for idx, mfrc in zip(idx_father[good_father],
                              m_frac_prg[good_father])]

                    idx = idx_father[np.argmax(macc_father)]
                    # Criterion 3
                    if abs(np.log10(atree[i-1]["m"]/t[idx]["m"])) > np.log10(max_dM_frac):
                        print("Sudden change in mass!")
                        idx=-2
                except:
                    idx=-2


            if idx == -2:
                # is it really a new galaxy??
                early_enough = nstep < nstep_early
                small_enough = atree[i-1]["m"] < m_small
                if early_enough or small_enough:
                    break
                else:
                    # No, it is a broken tree.
                    for istep_back in range(nstep_back_max):
                        # vicinity?
                        #xyz_guessed=[[]]*3
                        # i should NOT change
                        #for pind in range(3):
                        #    z = np.polyfit(atree["nstep"][max([0,i-6]):i-1],
                        #                atree["xp"][max([0,i-6]):i-1,pind], deg=2)
                        #    p = np.poly1d(z)
                        #    print(nstep_now)
                        #    xyz_guessed[pind] = p(nstep_now-istep_back)
                        xyz_guessed = atree["xp"][i-1,:]/self.aexps[nstep_now-1]*0.704
                        print("xyz guess", xyz_guessed)

                        # Use KDTree. -> Memory use?
                        # kdt_dist_upper could be calculated as v*dt.

                        #print(nstep_now-istep_back-2, kdt_dist_upper)
                        # kdt and ngals starts from 0.
                        neighbor_dist, neighbor_ind = self.kdts[nstep_now-istep_back-2].query(xyz_guessed,
                                                                       k=10,
                                                                       distance_upper_bound=kdt_dist_upper)
                        #print(neighbor_dist)
                        #print(np.isfinite(neighbor_dist))
                        neighbor_dist = neighbor_dist[np.isfinite(neighbor_dist)]
                        neighbor_ind = neighbor_ind[np.isfinite(neighbor_dist)]
                        #print(neighbor_dist, neighbor_ind)
                        if len(neighbor_dist) == 0:
                            # There is no hope..
                            continue
                        else:
                            print(np.sum(self.ngals[:nstep-istep_back])+neighbor_ind)
                            kdt_candidates = self.tree[np.sum(self.ngals[:nstep-istep_back])+neighbor_ind]
                            m_ratio = kdt_candidates["m"]/atree[i-1]["m"]
                            cvel_ratio = kdt_candidates["cvel"]/atree[i-1]["cvel"]
                        if len(neighbor_dist) == 1:
                            # final check
                            if (m_ratio) and ():
                                idx = kdt_candidates["idx"]
                            else:
                                # No hope
                                continue
                        elif len(neighbor_dist) > 1:
                            dist_norm = neighbor_dist/np.median(neighbor_dist)
                            # mass
                            # If it is missing, probably the mass is too small, not too large.
                            m_ratio_norm = m_ratio/np.median(m_ratio)
                            cvel_ratio_norm = cvel_ratio/np.median(cvel_ratio)
                            print("Scores:")
                            print(m_ratio, m_ratio_norm)
                            print(cvel_ratio, cvel_ratio_norm, dist_norm)
                            #ibest = np.argmin(m_ratio_norm+cvel_ratio_norm+dist_norm)
                            # Or,
                            #ibest = np.argmax(1./m_ratio_norm+1/.cvel_ratio_norm+1./dist_norm)
                            # Which is better?

                            idx = kdt_candidates["idx"][np.argmin(m_ratio_norm+cvel_ratio_norm+dist_norm)]
                            print("idx", idx)
                            break # get out of iteration.

            if idx < 1:
                # No prg FOR SURE!
                break

            if verbose:
                print("{}, M_son {:.2e}, M_now {:.2e}".format(idx, atree[i-1]["m"],t[idx]["m"]))

            atree[i]=t[idx]
            nouts.append(nstep)
        else:
            break

    return np.copy(atree[:i])


def extract_direct_full_tree(self, idx,
                             return_id=False,
                             return_macc=False,
                             max_dM_frac=5.0,
                             m_frac_min = 0.5,
                             verbose=False):
    """
    Extracts main progenitors from a TreeMaker tree.

    example
    -------
    >>> tt = tmtree.Tree("tree.dat")
    >>> atree = tt.extract_main_tree(12345)

    TODO
    ----
    It works, but the try - except clause is error-prone.
    Explicitly check the end of progenitor tree and make the function more predictable.

    """

    t = self.tree
    if return_id:
        fatherID = self.fatherID
    fatherIDx = self.fatherIDx
    fatherMass = self.fatherMass

    t_now = t[idx]
    nstep = t_now["nstep"]
    print(nstep)
    nouts = [nstep]
    atree = np.zeros(nstep + 1, dtype=t.dtype)
    atree[0] = t_now

    idx_prgs_alltime = [[idx]]

    if return_id:
        id_prgs_alltime = [[t[idx]["id"]]]
    if return_macc:
        macc_prgs_alltime = []
    for i in range(1, nstep + 1):
        if return_id:
            id_father  =  fatherID[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]# -1
        try:
            idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
            if len(idx_father) > 0:
                macc_father = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
                # In decending order of macc
                msort = np.argsort(macc_father)[::-1]
                idx_father = idx_father[msort]
                id_father = id_father[msort]
                macc_father = macc_father[msort]
                mass_father = np.array([t[idx]["m"] for idx in idx_father])
                m_frac_prg = atree[i-1]["m"] * (0.01*macc_father) / mass_father

                good_father = (m_frac_prg > m_frac_min) * (idx_father>0)

                idx_prgs_alltime.append(list(idx_father[good_father]))
                # list so to make it easy to remove later
                if return_id:
                    id_prgs_alltime.append(list(id_father[good_father]))

                if return_macc:
                    macc_prgs_alltime.append(macc_father[good_father])

                idx = idx_father[0]
                if idx < 1:
                    break
                if abs(np.log10(atree[i-1]["m"]/t[idx]["m"])) > np.log10(max_dM_frac):
                    if verbose:
                        print("{}, M_son {:.2e}, M_now {:.2e}".format(idx, atree[i-1]["m"],t[idx]["m"]))
                        print("Sudden change in mass!")
                    break
                t_father=t[idx]
                atree[i]=t_father
                nouts.append(nstep)
            else:
                break
        except:
            break
    if return_macc:
        if return_id:
            return atree, idx_prgs_alltime, id_prgs_alltime, macc_prgs_alltime
        else:
            return atree, idx_prgs_alltime, macc_prgs_alltime
    else:
        if return_id:
            return atree, idx_prgs_alltime, id_prgs_alltime
        else:
            return atree, idx_prgs_alltime


def plot_tree(axs, tree, i,j, alpha=0.3, sscale=1e-8, nnza=None, cmap="hsv"):
    colormap = cm.get_cmap(cmap)
    axs[0][0].scatter(tree["xp"][:,0],tree["xp"][:,1],
                      alpha=alpha,
                      s=tree["m"]*sscale,
                      #c=tree["idx"][0]%256,
                      #cmap=cmap,vmin=0, vmax=255,
                      label="{}-{}".format(i,j))
    axs[0][1].scatter(tree["xp"][:,1],tree["xp"][:,2],
                      s=tree["m"]*sscale,
                      #c=tree["idx"][0]%256,
                      #cmap=cmap,vmin=0, vmax=255,
                      alpha=alpha)
    axs[1][0].scatter(tree["xp"][:,2],tree["xp"][:,0],
                      s=tree["m"]*sscale,
                      #c=tree["idx"][0]%256,
                      #cmap=cmap,vmin=0, vmax=255,
                      alpha=alpha)
    if nnza is not None:
        xtime =nnza.step2nout(tree["nstep"])
    else:
        xtime = tree["nstep"]
    axs[1][1].plot(xtime, np.log10(tree["m"]), label="{}-{}".format(i,j),
                      color=colormap((tree["idx"][0]%256)/256.))
                      #vmin=0, vmax=255)
    #axs[1][1].scatter(xtime, np.log10(tree["m"]), s=5,
    #                  c=tree["idx"][0]%256,
    #                  cmap=cmap,
    #                  vmin=0, vmax=255)

def line_scatter(ax, x,y, c=None, cmap="hsv", s=5):
    ax.plot(x,y, c=c, vmin=0, vmax=255, cmap=cmap)
    ax.scatter(x,y,s=s, c=c, vmin=0, vmax=255, cmap=cmap)

def plot_tree_detail(axs, tree, i,j, alpha=0.5, nnza=None):
    if nnza is not None:
        xtime =nnza.step2nout(tree["nstep"])
    else:
        xtime = tree["nstep"]
    axs[0][0].scatter(tree["xp"][:,0],tree["xp"][:,1],
                      alpha=alpha,
                      s=tree["m"]*sscale,
                      label="{}-{}".format(i,j))
    axs[1][0].scatter(tree["xp"][:,1],tree["xp"][:,2],
                      s=tree["m"]*sscale,
                      alpha=alpha)
    axs[2][0].scatter(tree["xp"][:,2],tree["xp"][:,0],
                      s=tree["m"]*sscale,
                      alpha=alpha)
    line_scatter(axs[0][1],xtime,tree["vp"][:,0])
    line_scatter(axs[1][1],xtime,tree["vp"][:,1])
    line_scatter(axs[2][1],xtime,tree["vp"][:,2])
    line_scatter(axs[0][2],xtime,np.log10(np.abs(tree["ep"])))
    line_scatter(axs[1][2],xtime,np.log10(tree["ek"]))
    line_scatter(axs[2][2],xtime,np.log10(tree["et"]))
    line_scatter(axs[0][3],xtime,tree["spin"])
    line_scatter(axs[1][3],xtime,tree["cvel"])
    line_scatter(axs[2][3],xtime,np.log10(tree["m"]))


def check_tree(adp,
               save=True,
               nstep_min=0,
               detail=False,
               pos_diff=False,
               sscale=1e-8,
               nnza=None):
    """
        pos_diff is not working yet.
    """
    main = adp[0].pop(0)
    if pos_diff:
        sats = copy(adp)
        sats["x"]-=main["x"]
    else:
        sats = adp
    if detail:
        fig, axs = plt.subplots(3,4)
        fig.set_size_inches(12,8)
        plot_tree_detail(axs, main, 0,0,sscale=sscale, nnza=nnza)
        axs[0][0].set_xlabel(" X - Y ")
        axs[1][0].set_xlabel(" Y - Z ")
        axs[2][0].set_xlabel(" Z - X ")
        axs[0][1].set_xlabel(" vx ")
        axs[1][1].set_xlabel(" vy ")
        axs[2][1].set_xlabel(" vz ")
        axs[0][2].set_xlabel(" ep ")
        axs[1][2].set_xlabel(" ek ")
        axs[2][2].set_xlabel(" et ")
        axs[0][3].set_xlabel(" spin ")
        axs[1][3].set_xlabel(" cvel ")
        axs[2][3].set_xlabel(" m ")
    else:
        fig, axs = plt.subplots(2,2)
        fig.set_size_inches(8,6)
        plot_tree(axs, main, 0,0,sscale=sscale, nnza=nnza)
    for i, sats_this in enumerate(sats):
        for j, sat in enumerate(sats_this):
            if sat["nstep"][0] < nstep_min:
                break
            if detail:
                plot_tree_detail(axs,sat,i,j, sscale=sscale, nnza=nnza)
            else:
                plot_tree(axs,sat,i,j, sscale=sscale, nnza=nnza)
    axs[0][0].legend(markerscale=2.)
    plt.tight_layout()
    plt.suptitle("{}".format(main["idx"][0]))
    if save:
        plt.savefig("tree_check_{}.png".format(main["idx"][0]), dpi=300)
    else:
        plt.show()
    adp[0].append(main) # put it back.

def get_son(tt, idx):
    return tt.sonIDx[tt.tree[idx]["s_ind"]:tt.tree[idx]["s_ind"]+tt.tree[idx]["nsons"]]

def get_macc_father(tt, idx_org, merger_mass_frac_min=0.3, max_dM_frac=1e2):
    """

    """
    t = tt.tree
    fatherIDx = tt.fatherIDx
    fatherMass = tt.fatherMass
    if t[idx_org]["nprgs"] > 0:
        idx_father = fatherIDx[t["f_ind"][idx_org]:t["f_ind"][idx_org]+t["nprgs"][idx_org]]
        i_ok = idx_father > 0
        if sum(i_ok) > 0:
            idx_father = idx_father[i_ok]
            macc_father = fatherMass[t["f_ind"][idx_org]:t["f_ind"][idx_org]+t["nprgs"][idx_org]][i_ok]
            # In decending order of macc
            mass_father = np.array([t[fidx]["m"] for fidx in idx_father])
            m_frac_prg = t[idx_org]["m"] * (0.01*macc_father) / mass_father
            good_father = (m_frac_prg > merger_mass_frac_min) * (idx_father>0)
            if sum(good_father) == 0:
                return -1
            macc_father = macc_father[good_father]
            idx_father = idx_father[good_father]
            idx = idx_father[np.argmax(macc_father)]
            return idx

def macc_this_son(tt, idx_son, idx_father):
    i_father = np.where(idx_father == tt.fatherIDx[tt.tree["f_ind"][idx_son]:tt.tree["f_ind"][idx_son]+tt.tree["nprgs"][idx_son]])[0]
    # should be one.
    if len(i_father)>0:
        return tt.fatherMass[tt.tree["f_ind"][idx_son]+i_father][0], tt.tree[idx_son]["m"]

def filter_false_prg(tt,maintree, idx_prgs_alltime, mfrac_limit=50):
    for i, idx_prgs_now in enumerate(idx_prgs_alltime):
        if i==0:
            # i==0 = mainroot
            continue
        idx_son=maintree[i-1]["idx"]
        for idx in idx_prgs_now:
            son_macc, son_mass = macc_this_son(tt,idx_son, idx)
            # Fraction of father mass transferred to the son.
            if son_macc*son_mass/tt.tree[idx]["m"] < mfrac_limit:
                idx_prgs_now.remove(idx)
