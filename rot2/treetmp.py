import numpy as np
import matplotlib.pyplot as plt

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
                      max_dM_frac=5.0,
                      merger_mass_frac_min=0.5,
                      verbose=False):
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

    for i in range(1, nstep + 1):
        #print(i)
        idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
        i_ok = idx_father > 0
        if sum(i_ok) > 0:
            idx_father = idx_father[i_ok]
            macc_father = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]][i_ok]
            # In decending order of macc
            mass_father = np.array([t[idx]["m"] for idx in idx_father])
            m_frac_prg = atree[i-1]["m"] * (0.01*macc_father) / mass_father

            good_father = (m_frac_prg > merger_mass_frac_min) * (idx_father>0)
            if sum(good_father) == 0:
                break
            macc_father = macc_father[good_father]
            idx_father = idx_father[good_father]
            #idx = idx_father(np.argsort(macc_father)[::-1])
            idx = idx_father[np.argmax(macc_father)]# -1
            # Criterion 3
            if verbose:
                print(atree[i-1]["m"], t[idx]["m"])
            if abs(np.log10(atree[i-1]["m"]/t[idx]["m"])) > np.log10(max_dM_frac):
                print("{}, M_son {:.2e}, M_now {:.2e}".format(idx, atree[i-1]["m"],t[idx]["m"]))
                print("Sudden change in mass!")
                break

            if verbose:
                print("\n", t[idx]["m"],"step",i, "macc",mass_father)
                print("% sat", m_frac_prg)
                print("idx, m", [t[idx]["idx"] for idx in idx_father],[t[idx]["m"] for idx in idx_father])

            if idx < 1:
                break
            atree[i]=t[idx]
            nouts.append(nstep)
        else:
            break

    return np.copy(atree[:i])


def extract_direct_full_tree(self, idx,
                             return_id=False,
                             return_macc=False,
                             max_dM_frac=5.0,
                             merger_mass_frac_min = 0.5,
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

                good_father = (m_frac_prg > merger_mass_frac_min) * (idx_father>0)

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


def plot_tree(axs, tree, i,j, alpha=0.5, sscale=1e-8):
    axs[0][0].scatter(tree["xp"][:,0],tree["xp"][:,1],
                      marker='o',
                      #facecolors="none",
                      #edgecolors=np.random.random(),
                      alpha=alpha,
                      s=tree["m"]*sscale,
                      label="{}-{}".format(i,j))
    axs[0][1].scatter(tree["xp"][:,1],tree["xp"][:,2],
                      s=tree["m"]*sscale,
                      alpha=alpha)
    axs[1][0].scatter(tree["xp"][:,2],tree["xp"][:,0],
                      s=tree["m"]*sscale,
                      alpha=alpha)
    axs[1][1].plot(tree["nstep"], np.log10(tree["m"]), label="{}-{}".format(i,j))
    axs[1][1].scatter(tree["nstep"], np.log10(tree["m"]), s=5)


def plot_tree_detail(axs, tree, i,j, alpha=0.5, sscale=1e-8):
    axs[0][0].scatter(tree["xp"][:,0],tree["xp"][:,1],
                      alpha=alpha,
                      s=tree["m"]*sscale,
                      label="{}-{}".format(i,j))
    axs[0][1].scatter(tree["xp"][:,1],tree["xp"][:,2],
                      s=tree["m"]*sscale,
                      alpha=alpha)
    axs[0][2].scatter(tree["xp"][:,2],tree["xp"][:,0],
                      s=tree["m"]*sscale,
                      alpha=alpha)
    axs[1][0].plot(tree["nstep"],tree["vp"][:,1])
    axs[1][0].scatter(tree["nstep"],tree["vp"][:,1],s=5)
    axs[1][1].plot(tree["nstep"],tree["vp"][:,2])
    axs[1][1].scatter(tree["nstep"],tree["vp"][:,2],s=5)
    axs[1][2].plot(tree["nstep"],tree["vp"][:,0])
    axs[1][2].scatter(tree["nstep"],tree["vp"][:,0],s=5)
    axs[2][0].plot(tree["nstep"], np.log10(tree["spin"]))
    axs[2][1].plot(tree["nstep"], np.log10(tree["ek"]))
    axs[2][2].plot(tree["nstep"], np.log10(tree["m"]), label="{}-{}".format(i,j))
    axs[2][2].scatter(tree["nstep"], np.log10(tree["m"]), s=5)


def check_tree(adp, save=True, nstep_min=0, detail=False):
    main = adp[0].pop(0)
    sats = adp
    if detail:
        fig, axs = plt.subplots(3,3)
        fig.set_size_inches(12,8)
        plot_tree_detail(axs, main, 0,0)
        axs[0][0].set_xlabel(" X - Y ")
        axs[0][1].set_xlabel(" Y - Z ")
        axs[0][2].set_xlabel(" Z - X ")
        axs[1][0].set_xlabel(" vx ")
        axs[1][1].set_xlabel(" vy ")
        axs[1][2].set_xlabel(" vz ")
        axs[2][0].set_xlabel(" spin ")
        axs[2][1].set_xlabel(" ek ")
        axs[2][2].set_xlabel(" m ")
    else:
        fig, axs = plt.subplots(2,2)
        fig.set_size_inches(8,6)
        plot_tree(axs, main, 0,0)
    for i, sats_this in enumerate(sats):
        for j, sat in enumerate(sats_this):
            if sat["nstep"][0] < nstep_min:
                break
            if detail:
                plot_tree_detail(axs,sat,i,j, sscale=1e-8)
            else:
                plot_tree(axs,sat,i,j, sscale=1e-8)
    axs[0][0].legend(markerscale=2.)
    plt.tight_layout()
    plt.suptitle("{}".format(main["idx"][0]))
    if save:
        plt.savefig("tree_check_{}.png".format(main["idx"][0]), dpi=300)
    else:
        plt.show()

