import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import cm
from copy import copy

def build_kdtree(tt, H0):
    from scipy.spatial import cKDTree
    tt.kdts = []
    tt.idx_last =[]
    for i in range(tt.nsteps):
        td = tt.tree[tt.tree["nstep"]==i+1]
        tt.kdts.append(cKDTree(td["xp"]/tt.aexps[i]*(0.01*H0)))
        tt.idx_last.append(td[0]["idx"])
        #print(td[0]["idx"], td[-1]["idx"])
        #print(td[0]["id"], td[-1]["id"])
    tt.idx_last = np.array(tt.idx_last)


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
                    #print("get_all_trees", j, sat)
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


def extract_main_tree_try_fix(self, idx,
                      mmin=3.3e8,
                      max_dM_frac=50.0,
                      m_frac_min=0.5,
                      verbose=False,
                      kdt_dist_upper = 0.4):
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
    nstep_back_max = 3 # At most 5 steps backwards.
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
    if verbose:
        print("e----xtract_main_tree started at", atree[0]["nstep"])
    for i in range(1, nstep + 1):
        if istep_back > 1:
            istep_back -=1
            print("skipping ", i)
            # skip as many as istep_back
            continue
        connect = 0
        nstep_now = atree[i-1]["nstep"]
        idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
        i_ok = idx_father > 0
        if sum(i_ok) > 0:
            idx_father = idx_father[i_ok]
            macc_father = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]][i_ok]
            # In decending order of macc
            mass_father = np.array([t[fidx]["m"] for fidx in idx_father])
            m_frac_prg = atree[i-1]["m"] * (0.01*macc_father) / mass_father

            good_father = (m_frac_prg > m_frac_min)# * (idx_father>0)
            if sum(good_father) > 1:
                print("\n Father candidates before")
                [print("M_father_frac{:.2f}%  M_son_frac {:.2f}".format(100*mfrc, mac)) for mfrc, mac in zip(m_frac_prg, macc_father)]

            if sum(good_father) == 0:
                idx=-2
            else:
                #print("1  ", macc_father)
                macc_father = macc_father[good_father]
                #print("2   ",macc_father)
                idx_father = idx_father[good_father]

                if verbose:
                    print("\n Father candidates")
                    [print("{} {:.2f}%".format(idx, 100*mfrc)) for idx, mfrc in zip(idx_father,
                          m_frac_prg[good_father])]

                idx = idx_father[np.argmax(macc_father)]
                #print("iDX = ", idx)
                # Criterion 3
                if abs(np.log10(atree[i-1]["m"]/t[idx]["m"])) > np.log10(max_dM_frac):
                    print("Sudden change in mass!")
                    idx=-2
                #except:
                #    print("Except")
                #    idx=-2


            if idx == -2:
                print("IDX = -2  at nstep = ", nstep_now)
                # is it really a new galaxy??
                early_enough = nstep_now < nstep_early
                small_enough = atree[i-1]["m"] < m_small
                if early_enough or small_enough:
                    print("early or small", early_enough, small_enough)
                    break
                else:
                    # No, it is a broken tree.
                    # Do not attemp to find in the previous snapshot.
                    # TreeMaker confirmed that there is no hope.
                    if i < 3:
                        break
                    print("!!!!!!!!!", i)
                    ipoly=min([i,5])
                    print("ipoly", ipoly)
                    z = np.polyfit(atree["nstep"][i-ipoly:i-1],
                        atree["xp"][i-ipoly:i-1,0]/self.aexps[nstep_now-1:nstep_now-ipoly:-1]*0.704,
                        deg=2)
                    pxx = np.poly1d(z)
                    z = np.polyfit(atree["nstep"][i-ipoly:i-1],
                        atree["xp"][i-ipoly:i-1,1]/self.aexps[nstep_now-1:nstep_now-ipoly:-1]*0.704,
                        deg=2)
                    pyy = np.poly1d(z)
                    z = np.polyfit(atree["nstep"][i-ipoly:i-1],
                        atree["xp"][i-ipoly:i-1,2]/self.aexps[nstep_now-1:nstep_now-ipoly:-1]*0.704,
                        deg=2)
                    pzz = np.poly1d(z)



                    for istep_back in range(1,nstep_back_max+1):
                        # vicinity?
                        xyz_guessed = [pxx(nstep_now-istep_back),
                                       pyy(nstep_now-istep_back),
                                       pzz(nstep_now-istep_back)]
                        #xyz_guessed = atree["xp"][i-1,:]/self.aexps[nstep_now-1]*0.704
                        #xyz_guessed = [ -2.47599665, -48.05697326,  46.40619693]
                        #

                        # kdt_dist_upper could be calculated as v*dt.

                        #print(nstep_now-istep_back-2, kdt_dist_upper)
                        # kdt and ngals starts from 0.
                        neighbor_dist, neighbor_ind = self.kdts[nstep_now-istep_back-1].query(xyz_guessed,
                                                                       k=10*(1+istep_back),
                                                                       distance_upper_bound=kdt_dist_upper*(1+0.5*istep_back))
                        #answer = self.tree[self.idx_last[nstep_now-istep_back-1]-1+147266]
                        #print(answer["nstep"],
                        #      answer["idx"],
                        #      answer["xp"]/self.aexps[nstep_now-istep_back]*0.704)
                        #print(neighbor_dist)
                        #print(np.isfinite(neighbor_dist))
                        good = np.isfinite(neighbor_dist)
                        neighbor_dist = neighbor_dist[good]
                        neighbor_ind = neighbor_ind[good]
                        #print(neighbor_dist, neighbor_ind)
                        if len(neighbor_dist) == 0:
                            print("No neighbor found...")
                            # There is no hope..
                            continue
                        else:
                            #print(neighbor_ind)
                            #print(np.sum(self.ngals[:nstep_now-istep_back])+neighbor_ind)
                            kdt_candidates = self.tree[self.idx_last[nstep_now-istep_back-1]+neighbor_ind]
                            #print(kdt_candidates["nstep"][0])
                            #print(kdt_candidates["idx"])
                            #print(kdt_candidates["xp"])
                            m_ratio = kdt_candidates["m"]/atree[i-1]["m"]
                            cvel_ratio = np.abs(np.log2(kdt_candidates["cvel"]/atree[i-1]["cvel"]))
                            i_fine = np.where((m_ratio > 1/3) * (m_ratio < 5) * (cvel_ratio < 1))[0]
                            #kdt_candidates=kdt_candidates[i_fin]
                        if len(i_fine) == 1:
                            # final check
                            #if (m_ratio) and ():
                            idx = kdt_candidates["idx"][0]
                            #else:
                                # No hope
                            break
                        elif len(i_fine) > 1:
                            #dist_norm = neighbor_dist[i_fine]/np.mean(neighbor_dist[i_fine])
                            # mass
                            # If it is missing, probably the mass is too small, not too large.
                            m_ratio_std = (m_ratio[i_fine]-np.mean(m_ratio[i_fine]))/std(np.mean(m_ratio[i_fine]))
                            cvel_ratio_norm = cvel_ratio[i_fine]/np.mean(cvel_ratio[i_fine])
                            print("Scores:")
                            #print(m_ratio[i_fine])
                            print(m_ratio_std)
                            print(cvel_ratio[i_fine])
                            print(cvel_ratio_norm)
                            print(kdt_candidates[i_fine]["idx"])
                            scores = np.abs(np.log2(m_ratio_norm))+cvel_ratio_norm
                            ibest = np.argmin(scores)
                            #scores = 1./np.abs(np.log2(m_ratio_norm))+1./cvel_ratio_norm#+1./dist_norm
                            #ibest = np.argmax(scores)
                            print("final scores", scores)
                            #if scores[ibest] < 3.0 * np.max(scores[scores < scores[ibest]]):
                            print(m_ratio[ibest], cvel_ratio[ibest], neighbor_dist[ibest])
                            idx = kdt_candidates["idx"][i_fine[ibest]]
                            print("idx", idx)
                            if scores[ibest] > 1./2. * np.min(scores[scores > scores[ibest]]):
                                print("Not sure..")
                                continue
                            # Or,
                            #ibest = np.argmax(1./m_ratio_norm+1/.cvel_ratio_norm+1./dist_norm)
                            # Which is better?
                            connect=istep_back
                            break # get out of iteration.
            if idx < 1:
                # No prg FOR SURE!
                break

            if verbose:
                print("{}, M_son {:.2e}, M_now {:.2e}".format(idx, atree[i-1]["m"],t[idx]["m"]))
            if nstep_now == 179:
                print("idx", idx)
                print("istep_back", istep_back)
                print("i",i)
            if connect > 0:
                atree[i+connect-1]=t[idx]
            else:
                atree[i]=t[idx]
            nouts.append(nstep)
        else:
            break
    print("This tree is DONE at {}\n\n".format(nstep_now))

    return np.copy(atree[:i])

def extract_main_tree(self, idx,
                      mmin=3.3e8,
                      max_dM_frac=50.0,
                      m_frac_min=0.5,
                      verbose=False,
                      kdt_dist_upper = 0.4):
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
    nstep_back_max = 3 # At most 5 steps backwards.
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

    for i in range(1, nstep + 1):
        nstep_now = atree[i-1]["nstep"]
        idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
        i_ok = idx_father > 0
        if sum(i_ok) > 0:
            idx_father = idx_father[i_ok]
            macc_father = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]][i_ok]
            # In decending order of macc
            mass_father = np.array([t[fidx]["m"] for fidx in idx_father])
            m_frac_prg = atree[i-1]["m"] * (0.01*macc_father) / mass_father

            good_father = (m_frac_prg > m_frac_min)# * (idx_father>0)
            #if sum(good_father) > 1:
            #    print("\n Father candidates before")
            #    [print("M_father_frac{:.2f}%  M_son_frac {:.2f}".format(100*mfrc, mac)) for mfrc, mac in zip(m_frac_prg, macc_father)]

            if sum(good_father) == 0:
                idx=-2
            else:
                #print("1  ", macc_father)
                macc_father = macc_father[good_father]
                #print("2   ",macc_father)
                idx_father = idx_father[good_father]

                if verbose:
                    print("\n Father candidates")
                    [print("{} {:.2f}%".format(idx, 100*mfrc)) for idx, mfrc in zip(idx_father,
                          m_frac_prg[good_father])]

                idx = idx_father[np.argmax(macc_father)]
                #print("iDX = ", idx)
                # Criterion 3
                if abs(np.log10(atree[i-1]["m"]/t[idx]["m"])) > np.log10(max_dM_frac):
                    print("Sudden change in mass!")
                    idx=-2

            if idx < 1:
                # No prg FOR SURE!
                break

            if verbose:
                print("{}, M_son {:.2e}, M_now {:.2e}".format(idx, atree[i-1]["m"],t[idx]["m"]))

            atree[i]=t[idx]
            nouts.append(nstep)
        else:
            break
    #print("This tree is DONE at {}\n\n".format(nstep_now))

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

def line_scatter(ax, x,y, c=None, cmap="hsv", color=None, s=2):
    ax.plot(x,y, c=c)
    ax.scatter(x,y,s=s, c=c, vmin=0, vmax=255, cmap=cmap)

def plot_tree_detail(axs, tree, main=None,
                     i=0, j=0,
                     alpha=0.5, nnza=None, sscale=1e-8):
    if nnza is not None:
        xtime =nnza.step2nout(tree["nstep"])
    else:
        xtime = tree["nstep"]

    if main is not None:
        main_part = main[np.in1d(main["nstep"], tree["nstep"])]
        pos=main_part["xp"]-tree["xp"]
    else:
        pos = tree["xp"]

    axs[0][0].scatter(pos[:,0],pos[:,1],
                      alpha=alpha,
                      s=tree["m"]*sscale,
                      label="{}-{}".format(i,j))
    axs[1][0].scatter(pos[:,1],pos[:,2],
                      s=tree["m"]*sscale,
                      alpha=alpha)
    axs[2][0].scatter(pos[:,2],pos[:,0],
                      s=tree["m"]*sscale,
                      alpha=alpha)
    line_scatter(axs[0][1],xtime,tree["vp"][:,0])
    line_scatter(axs[1][1],xtime,tree["vp"][:,1])
    line_scatter(axs[2][1],xtime,tree["vp"][:,2])
    line_scatter(axs[0][2],xtime,tree["lp"][:,0])
    line_scatter(axs[1][2],xtime,tree["lp"][:,1])
    line_scatter(axs[2][2],xtime,tree["lp"][:,2])
    line_scatter(axs[0][3],xtime,np.abs(tree["rho_0"]))
    line_scatter(axs[1][3],xtime,tree["rs"]) # = Rs
    line_scatter(axs[2][3],xtime,np.log10(tree["ek"]/tree["m"]))
    line_scatter(axs[0][4],xtime,tree["spin"])
    line_scatter(axs[1][4],xtime,tree["cvel"])
    line_scatter(axs[2][4],xtime,np.log10(tree["m"]))


def check_tree(adp,
               out_dir="./",
               save=True,
               suffix="org",
               nstep_min=0,
               nstep_max=1e5,
               figure_type="regular",
               pos_diff=False,
               sscale=1e-8,
               nnza=None,
               cmap="hsv"):
    """
        pos_diff is not working yet.
    """
    main = adp[0].pop(0)
    if pos_diff:
        sats = copy(adp)
        sats["x"]-=main["x"]
    else:
        sats = adp
    if figure_type=="detailed":
        fig, axs = plt.subplots(3,5)
        fig.set_size_inches(20,12)
        if nstep_max > main["nstep"][0]:
            plot_tree_detail(axs, main, i=0,j=0,sscale=sscale, nnza=nnza)

        axs[0][0].set_xlabel(" X - Y ")
        axs[1][0].set_xlabel(" Y - Z ")
        axs[2][0].set_xlabel(" Z - X ")
        axs[0][1].set_xlabel(" vx ")
        axs[1][1].set_xlabel(" vy ")
        axs[2][1].set_xlabel(" vz ")
        axs[0][2].set_xlabel(" lx ")
        axs[1][2].set_xlabel(" ly ")
        axs[2][2].set_xlabel(" lz ")
        axs[0][3].set_xlabel(" rho_0 ")
        axs[1][3].set_xlabel(" Rs ")
        axs[2][3].set_xlabel(" ek ")
        axs[0][4].set_xlabel(" spin ")
        axs[1][4].set_xlabel(" cvel ")
        axs[2][4].set_xlabel(" m ")
    elif figure_type=="regular":
        fig, axs = plt.subplots(2,2)
        fig.set_size_inches(8,6)
        plot_tree(axs, main, i=0,j=0, sscale=sscale, nnza=nnza)
    elif figure_type=="simple":
        colormap = cm.get_cmap(cmap)
        fig, axs = plt.subplots()
        fig.set_size_inches(8,6)
        xtime = main["nstep"]
        axs.plot(xtime, np.log10(main["m"]), label="{}-{}".format(0,0),
                      color=colormap((main["idx"][0]%256)/256.))

    for i, sats_this in enumerate(sats):
        for j, sat in enumerate(sats_this):
            if sat["nstep"][0] < nstep_min or sat["nstep"][0] > nstep_max:
                break
            if figure_type=="detailed":
                plot_tree_detail(axs,sat,main=main,i=i,j=j, sscale=sscale, nnza=nnza)
            elif figure_type=="regular":
                plot_tree(axs,sat,i,j, sscale=sscale, nnza=nnza)
            elif figure_type=="simple":
                xtime = sat["nstep"]
                axs.plot(xtime, np.log10(sat["m"]), label="{}-{}".format(i,j),
                              color=colormap((sat["idx"][0]%256)/256.))

    if figure_type!="simple": axs[0][0].legend(markerscale=2.)
    plt.tight_layout()
    plt.suptitle("{}".format(main["idx"][0]))
    if save:
        plt.savefig(out_dir+"tree_check_{}_{}_{}.png".format(main["idx"][0], suffix, figure_type), dpi=300)
    else:
        plt.show()
    plt.close()
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


def singletree2txt(f, tree,fields_save,fields_fmt):
    data=[]
    fmts=[]
    for field, fmt in zip(fields_save, fields_fmt):
        if tree[field].ndim == 2:
            for i in range(3):
                data.append(tree[field][:,i])
                fmts.append(fmt)
        else:
            data.append(tree[field])
            fmts.append(fmt)
    #print(np.c_[data].shape)
    np.savetxt(f, np.column_stack(data), fmt=fmts, footer="\n")


def save_adp(adp, out_dir="./", suffix=""):
    import os
    if not os.path.isdir(out_dir):os.mkdir(out_dir)
    fields_save=["nstep", "id","idx","m","xp","vp","lp","mvir","rvir","tvir","cvel","rho_0","rs","ek","ep","et"]
    fields_fmt =["%d", "%d", "%d","%.4e","%.5f","%.5f","%.5f","%.4e","%.5f","%.5f","%.5f","%.5f","%.5f","%.5f","%.5f","%.5f"]
    f = open(out_dir+"{}{}.txt".format(adp[0][0]["idx"][0],suffix), "ab")
    header='      '.join(fields_save) + "\n"
    #np.savetxt(f, np.array([1]), header=header.encode('utf-8'))
    f.write(header.encode('utf-8'))
    for this_sats in adp:
        for sat in this_sats:
            singletree2txt(f, sat, fields_save, fields_fmt)
            f.write("\n".encode('utf-8'))
        f.write("\n".encode('utf-8'))
    f.close()
