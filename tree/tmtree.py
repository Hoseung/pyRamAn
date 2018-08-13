# -*- coding: utf-8 -*-

"""
From tree.fits file, select halos by mass cut.
Created on Wed Jan 21 11:45:42 2015

@author: hoseung
"""
import numpy as np
from utils.io_module import read_fortran, skip_fortran
from utils.hagn import Nnza
from utils import cosmology
import os
from os.path import isfile
from tree import halomodule

class Tree():
    def __init__(self, fn=None,
                     wdir='./',
                     load=True,
                     nout_now=None,
                     info=None,
                     BIG_RUN=True,
                     is_gal=False,
                     load_info=True):
        """
        NH tree => BIG_RUN == False
        """
        self.fn = None
        self.tree = None
        self.fatherID = None
        self.fatherIDx = None
        self.fatherMass = None
        self.sonIDx = None
        self.n_all_halos = 0
        self.n_all_fathers = 0
        self.n_all_sons = 0
        self.is_gal = is_gal
        self.wdir = wdir
        self.aexps = None
        self.omega_ts = None
        self.age_univs = None
        self.nsteps = None
        self.info = info
        self.verbose=False

        if self.info is None and load_info:
            self.load_info()

        if fn is not None and isfile(fn):
            self.set_fn(fn)
        else:
            self.get_fn()

        if load:
            self.load(nout_now=nout_now, BIG_RUN=BIG_RUN)
        # Load nnza()
        try:
            self.nnza = Nnza(fname=self.fn.split("tree.dat")[0]+"nout_nstep_zred_aexp.txt")
        except:
            try:
                self.cal_nnza()
                self.nnza = Nnza(fname=self.fn.split("tree.dat")[0]+"nout_nstep_zred_aexp.txt")
            except:
                print("[warning] Can not load nnza")

    def load_info(self):
        from load.info import Info
        from utils import util

        nout_fi = util.get_last_snapshot(self.wdir)
        self.info = Info(base=self.wdir,nout=nout_fi)

    def cal_time(self, info=None):
        """
            Calculate look back time and update nnza["lbt"].
            Assuems now = the end of tree.
            If the last step is at z=1, then lbt=0Gyr and z=0.
        """
        if not hasattr(self, "nnza"):
            print("No nnza attribute found, Can't proceed")
            return

        if info is None:
           info = self.info

        if info is None:
            self.load_info()

        # Add error handling. Give useful information on error
        tc = cosmology.Timeconvert(info)
        z_now = min(self.nnza.nnza["zred"])
        self.nnza.nnza["lbt"] = tc.zred2gyr(self.nnza.nnza["zred"], z_now=z_now)
        print("Look back time calculation done")
        print("Current time : z={:.2f}".format(z_now))


    def dump(self, suffix="", force=False, protocol=-1):
        """
            Save tree data and class instance.
            After saving it, the data is gone!

            Parameters
            ----------
            force : default = False
                do not check before erasing data
            protocol : by default -1
                pickle protocol. Negative value = highest possible.
        """
        import pickle
        self.dump_files={}
        dump_temp=[]
        # Check if all set.
        for attr_name in ["tree", "fatherID", "fatherIDx", "fatherMass", "sonIDx"]:
            attr = getattr(self, attr_name)
            if attr is None:
                print("Nothing to save ", attr_name)
                return
            dump_temp.append(attr)

        fn = self.wdir + suffix + "data.npy"
        np.save(fn, dump_temp)
        self.dump_files.update({"data":fn})

        #if input("self.tree, self.father* data will be deleted, OK? [y/n]\n") == "y" or force:
        # empty the data for the moment.
        self.tree = None
        self.fatherID = None
        self.fatherIDx = None
        self.fatherMass = None
        self.sonIDx = None
        fn = self.wdir + suffix + "tree_meta_" + ["hal","gal"][self.is_gal]
        if not fn.endswith(".pickle"):
            fn += ".pickle"
        self.dump_files.update({"meta":fn})
        pickle.dump(self, open(fn, "wb"), protocol=protocol)
        # restore data
        self.tree, self.fatherID, self.fatherIDx, self.fatherMass, self.sonIDx = dump_temp

    def load_np(self, suffix="", protocol=-1):
        """
        Assuming the current (loaded) instance is the matching meta data.
        """
        Overwrite_ok = False
        for attr_name in ["tree", "fatherID", "fatherIDx", "fatherMass"]:
            attr = getattr(self, attr_name)
            if attr is not None and not Overwrite_ok:
                if input("Overwrite all data? [y/n]\n".format(attr_name)) != "y":
                    Overwrite_ok=True
        #            return
        self.tree, self.fatherID, self.fatherIDx, self.fatherMass, self.sonIDx = np.load(self.dump_files["data"])

    def get_fn(self):
        if self.is_gal:
            dir_mid = "GalaxyMaker/gal/"
        if not self.is_gal:
            dir_mid = "halo/DM/"
        fn = self.wdir + dir_mid + "tree.dat"

        if isfile(fn):
            self.fn = fn
        elif isfile(self.wdir + dir_mid.split("/")[0] + "/" + "tree.dat"):
            self.fn = self.wdir + dir_mid.split("/")[0] + "/" + "tree.dat"
            print("tree.dat found in another location", self.fn)
        else:
            #self.set_fn(None)
            print(fn, "is not found")


    def load(self, BIG_RUN=True, nout_now = None):
        """
            Parameters
            ----------

            nout_now (Deprecated) : int
                correct nstep to match with simulation nout.
                The first treebrick in the tree is usually between
                a few to a few tens snapshot.
                **
                nout - nstep relation is no more that simple due to missing nouts when building trees.
                    Thus it is better to keep a nout-nstep-zred-aexp(-lbt) matching list, which is hagn.Nnza().

            BIG_RUN : logical
                for large simulation, particle mass is excluded.
                .. But why 'np' of tree is excluded too?
                it's only one more column to ~40 coluns.
        """
        if self.fn is None:
            print("No tree file name is given")
            return
        from tree import cnt_tree

        print("BIG_RUN", BIG_RUN)
        self.n_all_halos, self.n_all_fathers, self.n_all_sons, self.nsteps = cnt_tree.count_tree(self.fn, int(BIG_RUN))

        if self.verbose:
            print(self.n_all_halos, self.n_all_fathers, self.n_all_sons, self.nsteps)
        self.fatherID, self.fatherIDx, self.sonIDx, \
        self.fatherMass, i_arr, f_arr, \
        self.aexps, self.omega_ts, self.age_univs = \
                cnt_tree.load_tree(self.fn, self.n_all_halos, \
                self.n_all_fathers, self.n_all_sons, int(BIG_RUN), self.nsteps)

        self.fatherIDx -=1
        # zred is omitted to reduce memory usage
        from load import dtypes
        dtype_tree = dtypes.get_tree_dtypes(BIG_RUN=BIG_RUN)

        tt = np.recarray(self.n_all_halos +1, dtype = dtype_tree)
        self.tree = tt

        tt["m"][1:] = f_arr[:,0] * 1e11
        tt["macc"][1:] = f_arr[:,1]
        tt["xp"][1:] = f_arr[:,2:5]
        tt["vp"][1:] = f_arr[:,5:8]
        tt["lp"][1:] = f_arr[:,8:11]
        tt["abc"][1:] = f_arr[:,11:15]
        tt["ek"][1:] = f_arr[:,15]
        tt["ep"][1:] = f_arr[:,16]
        tt["et"][1:] = f_arr[:,17]
        tt["spin"][1:] = f_arr[:,18]

        tt["rvir"][1:] = f_arr[:,19]
        tt["mvir"][1:] = f_arr[:,20]* 1e11
        tt["tvir"][1:] = f_arr[:,21]
        tt["cvel"][1:] = f_arr[:,22]
        tt["rho_0"][1:] = f_arr[:,23]
        tt["rs"][1:] = f_arr[:,24]

        tt["idx"][1:] = i_arr[:,0]
        tt["id"][1:] = i_arr[:,1]
        #tt["bushID"] = i_arr[:,2]
        #tt["st"] = i_arr[:,3]
        tt["level"][1:] = i_arr[:,4]
        tt["hosthalo"][1:] = i_arr[:,5]
        tt["hostsub"][1:] = i_arr[:,6]
        tt["nsub"][1:] = i_arr[:,7]
        tt["nextsub"][1:] = i_arr[:,8]
        tt["nprgs"][1:] = i_arr[:,9]
        if nout_now is not None:
            tt["nstep"][1:] = i_arr[:,10]
            #tt["nstep"][1:] += (nout_now - max(tt["nstep"]))
        else:
            tt["nstep"][1:] = i_arr[:,10]
        if not BIG_RUN:
            tt["np"][1:] = i_arr[:,11]
        tt["f_ind"][1:] = i_arr[:,12] -1 #
        tt["nsons"][1:] = i_arr[:,13] #
        tt["s_ind"][1:] = i_arr[:,14] -1 #


    def get_best_matched_desc(pids_now, gids, nout_next):
        gcat_with_pids = halomodule.Halo(nout=nout_next, return_id=gids, is_gal=True)
        n_matched=[]
        for i, pids in enumerate(gcat_with_pids.idlists):
            # pids = most_bound_partcles(pids)
            n_machted.append(np.len(np.intersect1d(pids, pids_now)))

        return np.argmax(n_matched), gcat_with_pids.idlists[np.argmax(n_matcehd)]


    def extract_main_tree_reverse(self, idx):
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
        fatherIDx = self.fatherIDx
        fatherMass = self.fatherMass

        t_now = t[idx]
        nstep = t_now["nstep"]
        nout_now = self.nnza.step2out(int(nstep))
        atree = np.zeros(nstep + 1, dtype=t.dtype)

        atree[0]=t_now
        nstep_max = t[-1]["nstep"]
        if nstep == nstep_max:
            return

        pids_now = halomodule.Halo(nout=nout_now, return_id=t_now, is_gal=True)[0]
        for i in range(1, nstep_max - nstep):
            try:
                nstep = t_now["nstep"]
                t_next = t[t["nstep"] == nstep+1]
                f_ind_first, f_ind_last = t_next["f_ind"][0], t_next["f_ind"][-1]+t_next["nprgs"]
                ind_matched = np.where(fatherIDx[f_ind_first:f_ind_last] == idx)[0]
                ind_next_step = np.searchsorted(np.cumsum(t_next["nprgs"]), ind_matched)
                all_desc = t_next[ind_next_step]
                i_best, pids_now = get_best_matched_desc(pids_now, all_desc["id"], nout_next = self.nnza.step2out(nstep+1))
                idx = all_desc["idx"][i_best]
                t_now = t[idx]

            except:
                break
            atree[i]=t_now

        return np.copy(atree[:i][::-1])


    def extract_direct_full_tree(self, idx, return_id=False):
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
            for i in range(1, nstep + 1):
                id_father  =  fatherID[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]# -1
                try:
                    idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]# -1
                    if len(idx_father) > 0:
                        idx_prgs_alltime.append(list(idx_father[idx_father>0]))
                        id_prgs_alltime.append(list(id_father[id_father>0]))
                        mass_father_transfer = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
                        mass_father = t[idx_father[idx_father>0]]["m"]
                        idx = idx_father[np.argmax(mass_father_transfer/mass_father)]
                        if idx < 1:
                            break
                        t_father=t[idx]
                        atree[i]=t_father
                        nouts.append(nstep)
                    else:
                        break
                except:
                    break
            return atree, idx_prgs_alltime, id_prgs_alltime
        else:
            for i in range(1, nstep + 1):
                try:
                    idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]# -1
                    if len(idx_father) > 0:
                        idx_prgs_alltime.append(list(idx_father[idx_father>0]))
                        mass_father_transfer = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
                        mass_father = t[idx_father[idx_father>0]]["m"]
                        idx = idx_father[np.argmax(mass_father_transfer/mass_father * (np.abs(np.log10(mass_father/t[idx]["m"])) < 0.2).astype(int))]
                        if idx < 1:
                            break
                        t_father=t[idx]
                        atree[i]=t_father
                        nouts.append(nstep)
                    else:
                        break
                except:
                    break

            return atree, idx_prgs_alltime

    def get_all_trees(self, idx_prgs_alltime,
					 skip_main=True,
					 filter_dup=True):
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


        """
        all_main_prgs=[]
        # loop over all nstep

        for j, satellite_roots in enumerate(idx_prgs_alltime):
            mainprgs=[]
            # loop over all satellites at each step
            for i,sat in enumerate(satellite_roots):
                if not skip_main or i!=0:
                    #print("sat ind", i)
                    mainprgs.append(self.extract_main_tree(sat))
            all_main_prgs.append(mainprgs)
            all_idxs_filter = []
            if filter_dup:
                if len(mainprgs) > 0:
                    for aa in mainprgs:
                        try:
                            all_idxs_filter.extend(aa["idx"][1:])
                        except:
                            print(mainprgs)
                            return mainprgs
                        # last idx MUST remain in the prgs.
                    for idxs in idx_prgs_alltime[j+1:]:
                        if len(idxs) > 1:
                            for idx in idxs[1:]:
                                if idx in all_idxs_filter:
                                    idxs.remove(idx)

        return all_main_prgs

    def extract_main_tree(self, idx):
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
        fatherIDx = self.fatherIDx
        fatherMass = self.fatherMass

        t_now = t[idx]
        nstep = t_now["nstep"]
        nouts = [nstep]
        atree = np.zeros(nstep + 1, dtype=t.dtype)
        atree[0] = t_now

        if nstep <= 1:
            return

        for i in range(1, nstep + 1):
            try:
                idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
                if len(idx_father) > 0:
                    mass_father = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
                    idx = idx_father[np.argmax(mass_father)]# -1
                    if idx < 1:
                        break
                    t_father=t[idx]
                    #idx = t_father["idx"]
                    atree[i]=t_father
                    nouts.append(nstep)
                else:
                    break
            except:
                break

        return np.copy(atree[:i])

    def cal_nnza(self):
        """

        ToDo: change name as get_nnza()
        """
        fdir = self.fn.split("tree.dat")[0]

        fsave = fdir+"nout_nstep_zred_aexp.txt"
        if os.path.isfile(fsave):
            print("File {} exists".format(fsave))
            #return
        f_tree_input = fdir + "input_TreeMaker.dat"
        if not os.path.isfile(f_tree_input):
            print("[warning] {} doest not exist".format(f_tree_input))
            #return

        with open(f_tree_input, "r") as f:
            nsteps = int(f.readline().split()[0])
            nouts = []
            for i,line in enumerate(f.readlines()):
                if "tree_bricks" in line:
                    nouts.append(int(line.split("tree_bricks")[1].split("'")[0]))

        nsteps = np.unique(self.tree["nstep"])[::-1] # decending order
        # remove nstep = 0 in an empty tree
        nsteps = nsteps[nsteps > 0]
        # Early snapshots in input_TreeMkaer.dat can be ignored.
        nnsteps = len(nsteps)
        nouts=np.array(nouts)[nnsteps::-1] # decending order
        aexps = self.aexps[nnsteps::-1] # decending order
        zreds = 1./aexps - 1
        np.savetxt(fsave, np.c_[nouts, nsteps, zreds, aexps], fmt=['%d', '%d', '%.9f', '%.9f'])

def load_fits(base=None, work_dir=None, filename="halo/TMtree.fits"):
    """
        Load tree.fits file and get the data table.
        use base always. Work_dir is for backward compatibility.
    """
    from astropy.io import fits
    from astropy.table import Table
    if base == None:
        base = work_dir

    if base == None:
        raise ValueError ("Please provide the parameter base=")

    data = fits.getdata(base + filename, 1)
    return Table(data)


def fix_nout(tt, nout_ini, nout_fi):
    nout_min_org = tt['NOUT'].min()
    nout_max_org = tt['NOUT'].max()

    nnout_org = nout_max_org - nout_min_org + 1
    nnout_new = nout_fi - nout_ini + 1

    assert (nnout_org == nnout_new), "number of nouts does not match"

    i_z_max = np.where(tt["Z"] == tt["Z"].max())[0]
    assert (tt["NOUT"][np.where(tt["Z"] == tt["Z"].max())[0]][0] == nout_max_org), "The snapshot number of highest redshift snapshot is not 0. Maybe you've already fixed nouts?"

    # OK. It's safe.
    tt["NOUT"] = nout_fi - tt["NOUT"]


def check_tree_complete(tree, nout_fi, nout_ini, halo_list):
    import numpy as np
    '''
    returns a list of halo IDs at nout_fi that with trees fully linked
    from nout_ini to nout_fi.

    '''
    # Make sure all parameters are given
    complete_list=np.zeros(len(halo_list), dtype=bool)
    #for ihal,halo in enumerate(halo_list):
    idxlist, idlist=get_main_prg(tree, halo_list, nout_fi=0, nout_ini=nout_ini)
    for i in range(len(halo_list)):
        if len(np.where(idxlist[i] == 0)[0]) > 1:
            complete = False
        elif len(np.where(idxlist[i] == 0)[0]) == 0:
            complete = True
        complete_list[i] = complete

    return complete_list, idlist[complete_list]




def check_tree_fig(tt, idx):
    import matplotlib.pyplot as plt
    aexp = tt.nnza.nnza["aexp"]
    atree = extract_direct_full_tree(tt,idx)[0]
    fig, axs = plt.subplots(2,2)
    axs=axs.ravel()

    i_ok = np.where(atree["m"] > 0)[0]
    axs[0].scatter(atree["xp"][i_ok,0]/aexp[i_ok], atree["xp"][i_ok,1]/aexp[i_ok], c=atree["nstep"][i_ok])
    axs[1].scatter(atree["nstep"][i_ok],np.log10(atree["m"][i_ok]), c=atree["nstep"][i_ok])
    axs[2].scatter(atree["vp"][i_ok,0], atree["vp"][i_ok,1], c=atree["nstep"][i_ok])

    axs[0].set_xlabel("pos x")
    axs[0].set_ylabel("pos y")
    axs[1].set_xlabel("nstep")
    axs[1].set_ylabel("log(m*)")
    axs[2].set_xlabel("vel x")
    axs[2].set_ylabel("vel y")

    plt.savefig("tree_{}.png".format(idx))
    plt.close()

def extract_direct_full_tree(self, idx, return_id=False):
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
        for i in range(1, nstep + 1):
            id_father  =  fatherID[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]# -1
            try:
                idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]# -1
                if len(idx_father) > 0:
                    idx_prgs_alltime.append(list(idx_father[idx_father>0]))
                    id_prgs_alltime.append(list(id_father[id_father>0]))
                    mass_father_transfer = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]
                    mass_father = t[idx_father[idx_father>0]]["m"]
                    idx = idx_father[np.argmax(mass_father_transfer/mass_father)]
                    if idx < 1:
                        break
                    t_father=t[idx]
                    atree[i]=t_father
                    nouts.append(nstep)
                else:
                    break
            except:
                break
        return atree, idx_prgs_alltime, id_prgs_alltime
    else:
        m_now = t["m"][idx]
        for i in range(1, nstep + 1):
            print("\n step {}".format(i))
            print("m_now = ", np.log10(m_now))
            if True:
                idx_father = fatherIDx[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]]# -1
                if len(idx_father) > 0:
                    idx_prgs_alltime.append(list(idx_father[idx_father>0]))
                    idx_ok = idx_father>0
                    mass_father_transfer = fatherMass[t["f_ind"][idx]:t["f_ind"][idx]+t["nprgs"][idx]][idx_ok]*m_now
                    print("mass_father_transfer", mass_father_transfer)
                    mass_father = t[idx_father[idx_ok]]["m"]
                    mass_ratio = mass_father_transfer/mass_father
                    #is_ok_mass = np.abs(np.log10(mass_father_transfer/t[idx]["m"])-2) < 0.2
                    is_ok_mass = np.abs(np.log10(mass_ratio)-2) < 0.2
                    print("mass_ratio", mass_ratio)
                    print("ok_mass", is_ok_mass)
                    print("product",mass_ratio * is_ok_mass)
                    if sum(is_ok_mass) == 0:
                        break
                    idx = idx_father[np.argmax(mass_ratio * is_ok_mass)]
                    print(np.argmax(mass_ratio * is_ok_mass), idx)
                    if idx < 1:
                        break
                    t_father=t[idx]
                    atree[i]=t_father
                    nouts.append(nstep)
                    m_now = 10**np.mean(np.log10(atree["m"][:i]))
                    print(atree["m"])
                else:
                    print("break 2")
                    break
            else:
                print("break 3")
                break

        return atree, idx_prgs_alltime
