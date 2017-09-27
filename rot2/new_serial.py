from utils import match as mtc
import rot2.new_serial_modules as nsm
import numpy as np
import pickle

dtype_finedata = [('nstep', '<i4'),
                  ('nout', '<i4'),
                  ('lbt', '<f8'),
                  ('idx', '<i4'),
                  ('id', '<i4'),
                  ('mstar', '<f8'),
                  ('pos', '<f8', (3,)),
                  ('vel', '<f8', (3,)),
                  ('rgal', '<f8'),
                  ('cvel', '<f8'),
                  ('P_tidal', '<f8'),
                  ("P_tidal_h", "<f8"),
                  ('d5', '<f8'),
                  ('d10', '<f8'),
                  ('d50', '<f8'),
                  ("mgas", "<f8"),
                  ("mgas_cold", "<f8"),
                  ("reff", "<f8"),
                  ("lambda_r", "<f8"),
                  ("eps", "<f8"),
                  ("vsig", "<f8"),
                  ("sfr01", "<f8"),
                  ("sfr05", "<f8"),
                  ("sfr1", "<f8"),
                  ("sfr_area", "<f8"),
                  ("lgas", "<f8", (3)),
                  ("lvec", "<f8", (3)),
                  ("nvec", "<f8", (3))]

dtype_merger_props = [("lbt", '<f8'),
                    ("nout", '<i8'),
                    ("id", '<i8'),
                    ("idx", '<i8'),
                    ("dist", "<f8"),
                    ("orbitang", '<f8'),
                    ("rel_pos", '<f8', (3,)),
                    ("rel_vel", '<f8', (3,)),
                    ("jorbit", '<f8', (3,)),
                    # Below are to be interpolated from measurements.
                    ("spinang", '<f8'),
                    ("spinmag", '<f8'),
                    ("mstar", '<f8'),
                    ("m_frac", '<f8'),
                    ("m_gas", '<f8'),
                    ("m_gas_cold", '<f8'),
                    ("rgal", '<f8'),
                    ("rgal_sum","<f8"),
                    ("reff_s","<f8")]

dtype_main_results = [("nstep", "<i4"),
                      ("nout", "<i4"),
                      ("lbt", "<f8"),
                      ("id", "<i4"),
                      ("idx", "<i4"),
                      ("mstar", "<f8"),
                      ("mgas", "<f8"),
                      ("mgas_cold", "<f8"),
                      ("rgal", "<f8"),
                      ("reff", "<f8"),
                      ("lambda_r", "<f8"),
                      ("eps", "<f8"),
                      ("vsig", "<f8"),
                      ("sfr01", "<f8"),
                      ("sfr05", "<f8"),
                      ("sfr1", "<f8"),
                      ("sfr_area", "<f8"),
                      ("lgas", "<f8", (3)),
                      ("pos", "<f8", (3)),
                      ("vel", "<f8", (3)),
                      ("lvec", "<f8", (3)),
                      ("nvec", "<f8", (3)),
                      ("P_tidal_h", "<f8")]

dtype_sat_results = [ ("nout", "<i4"),
                      ("lbt", "<f8"),
                      ("mstar", "<f8"),
                      ("mgas", "<f8"),
                      ("mgas_cold", "<f8"),
                      ("rgal", "<f8"),
                      ("reff", "<f8"),
                      ("lambda_r", "<f8"),
                      ("vsig", "<f8"),
                      ("sfr01", "<f8"),
                      ("sfr05", "<f8"),
                      ("sfr1", "<f8"),
                      ("sfr_area", "<f8"),
                      ("lgas", "<f8", (3)),
                      ("pos", "<f8", (3)),
                      ("vel", "<f8", (3)),
                      ("lvec", "<f8", (3)),
                      ("nvec", "<f8", (3))]

min_sattree_length = 21

class Maingal():
    def __init__(self, adp, nnza_all, nnouts, nnza_cell):
        self.i_last_main=0 # main_arr data input index
        self.fidx = adp[0][0]["idx"][0]
        # finedata is poppeds
        self.set_finedata(adp[0].pop(0), nnza_all)

        # Now Only sats are remaing.
        self.main_data=np.zeros(nnouts, dtype=dtype_main_results)
        self.main_data["nout"] = nnza_cell.nnza["nout"][:nnouts]
        self.main_data["lbt"] = nnza_cell.nnza["lbt"][:nnouts]
        # What if main gal measurements are not available at early nouts?
        self.mergers=[]
        self.sat_data=[]
        self.add_mergers_tree(adp, nnza_cell, nnouts)
        self.n_minor = 0
        self.n_major = 0
        self.nout_data = nnza_cell.nnza["nout"][:nnouts]

    def set_finedata(self, tree, nnza_all):
        tree = tree[np.where(tree["nstep"] > 0)[0]]
        self.finedata = np.zeros(len(tree["nstep"]), dtype=dtype_finedata)
        self.finedata["pos"] = tree["xp"]
        self.finedata["vel"] = tree["vp"]
        self.finedata["id"] = tree["id"]
        self.finedata["idx"] = tree["idx"]
        self.finedata["nstep"] = tree["nstep"]
        self.finedata["mstar"] = tree["m"] # to be
        self.finedata["lbt"] = nnza_all.a2b(self.finedata["nstep"],"nstep","lbt")
        self.finedata["nout"] = nnza_all.a2b(self.finedata["nstep"],"nstep","nout")


    def add_mergers_tree(self, adp, nnza_cell, nnouts):
        """
            Add merger_arr from sat fine merger trees.
        """
        # sat should be no longer than the main.
        #sattree = sattree[sattree["nstep"] > self.finedata["nstep"].min()]
        #mm = Merger(self.main_arr[self.main_arr["nout"] > min(sat_nouts)], results, sattree, self)
        # main longer than sat is not needed.
        # change field name.
        # nextsub -> nout
        for adps_now in adp:
            for sattree in adps_now:
                marr = cal_merger_props_tree(self, sattree)
                if marr is not None:
                    self.mergers.append(marr)

        # sat_data to hold sat measurements
        for mm in self.mergers:
            matched_nouts =np.intersect1d(mm["nout"],nnza_cell.nnza["nout"][:nnouts])
            sd = np.full(len(matched_nouts), np.nan, dtype=dtype_sat_results)
            sd["nout"]=matched_nouts[::-1]
            sd["lbt"]=nnza_cell.a2b(sd["nout"], "nout", "lbt")
            self.sat_data.append(sd)


def cal_merger_props_tree(this_gal, sattree,
                           mass_ratio_cut=1/50, min_len=15):
    """
        Calculate merger related quantities, this time with only tree information.

        Parameters
        ----------
        mass_ratio_cut: 1/50
            If the maximum mass ratio is below mass_ratio_cut, ignore.
        min_len : 15
            Sat tree shorter than 15 is ignored.
            (15 in fine tree ~= 1 in coarse tree)

        NOTE
        ----
            More through mass ratio test must be done onnce
            the beginning of merger is determined.
    """
    if len(sattree) < 15:
        return None

    # only overlapping part of trees.
    finedata = this_gal.finedata
    common_nsteps = np.intersect1d(finedata["nstep"], sattree["nstep"])[::-1]
    if len(common_nsteps) < min_sattree_length:
        return None

    main_tree_part = finedata[mtc.match_list_ind(finedata["nstep"],common_nsteps)]
    sattree = sattree[mtc.match_list_ind(sattree["nstep"],common_nsteps)]
    assert (np.all(main_tree_part["nstep"] == sattree["nstep"])), "main and tree are incompatible"

    # Even the trees are long and overlapping, measured part can only minimally overlap.
    # Abort in such cases.
    if this_gal.main_data["nout"][this_gal.main_data["nout"] > 0][-3] > main_tree_part["nout"].max():
        return None

    # If the sat is too small, throw away.
    if np.max(sattree["m"]/main_tree_part["mstar"]) < mass_ratio_cut:
        return None

    rel_pos = main_tree_part["pos"] - sattree["xp"]
    # if more than half is confused with the main tree,
    # throw it away.
    if np.sum(rel_pos[:,0] == 0) > len(rel_pos)/2:
        return None

    rel_vel = main_tree_part["vel"] - sattree["vp"]

    merger_arr = np.zeros(len(sattree), dtype=dtype_merger_props)
    merger_arr["nout"]=main_tree_part["nout"]
    # spin alignment
    #merger_arr["nstep"]=sattree["nstep"]

    merger_arr["id"]=sattree["id"]
    merger_arr["mstar"] = sattree["m"] # Just for checking
    merger_arr["idx"]=sattree["idx"]
    merger_arr["rel_pos"]=rel_pos
    merger_arr["rel_vel"]=rel_vel
    return merger_arr

def cal_merger_props_measurments(main_arr, sat_arr):
    merger_arr["orbitang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_tree_part["nvec"], j_orbital))
    merger_arr["reff_p"]=main_tree_part["reff"]
    merger_arr["reff_s"]=self.sat.reff
    merger_arr["spinmag"]=np.sqrt(np.einsum('...i,...i', self.sat.lvec,self.sat.lvec))
    merger_arr["spinang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_tree_part["nvec"], self.sat.nvec))
    return merger_arr



def add_main_result(self, ss):
    # 63 steps
    il=self.i_last_main
    self.measurements["nout"][il] = ss.nout
    self.measurements["nstep"][il] = ss.nstep
    self.measurements["id"][il] = ss.id
    self.measurements["idx"][il] = ss.idx
    self.measurements["mstar"][il] = ss.mstar
    self.measurements["reff"][il] = ss.reff
    self.measurements["pos"][il,0] = ss.xc
    self.measurements["pos"][il,1] = ss.yc
    self.measurements["pos"][il,2] = ss.zc
    self.measurements["vel"][il,0] = ss.vxc
    self.measurements["vel"][il,1] = ss.vyc
    self.measurements["vel"][il,2] = ss.vzc
    self.measurements["lambda_r"][il] = data.nout

    if hasattr(ss, "gas_results"):
        self.measurements["mgas"] = ss.gas_results["mgas_tot"]
        self.measurements["mgas_col"] = ss.gas_results["mgas_cold"]
        self.measurements["lgas"]= ss.gas_results["Ln_gas"]
    if ss.lvec is not None:
        # Todo
        # suppress nvec being nan when measuring it from stars.
        self.measurements["lvec"][il] = ss.lvec
        self.measurements["nvec"][il] = ss.nvec
        # If there is Lvec, there is rgal.
        self.measurements["rgal"][il] = ss.rgal
        try:
            self.measurements["lambda_r"][il]=ss.lambda_r[0]
        except:
            self.measurements["lambda_r"][il]=ss.lambda_r

    try:
        self.measurements["vsig"][il]=ss.vsig_results["V_sig"]
    except:
        pass

    if hasattr(ss, "sfr_results"):
        self.measurements["sfr01"][il]=ss.sfr_results["sfrs"][0]
        self.measurements["sfr05"][il]=ss.sfr_results["sfrs"][1]
        self.measurements["sfr1"][il]=ss.sfr_results["sfrs"][2]
        self.measurements["sfr_area"][il]=ss.sfr_results["area"]
#def finalize_main_result

from scipy.signal import savgol_filter
def sat_data_to_fine(this_gal, merger, sat_data,
                     nnza_all,
                     fields_to_interpol_from,
                     fields_to_interpol_to,
                     win_size_small=15,
                     win_size_large=21,
                     do_smooth=True):
    """
    Interpolate values needed to determining merger epoch only.
    Others are kept in sat_data and can be inerpolated when needed.

    Parameters
    ----------
     win_size_small=15
         smoothing window size for non-monotonic variables.
     win_size_large=21
         smoothing window size for monotonic variables.


    >>> sat_data.dtype
    ... [('nout', '<i4'), ('zred', '<f8'), ('mstar', '<f8'),
    ...  ('mgas', '<f8'), ('mgas_cold', '<f8'), ('rgal', '<f8'),
    ...  ('reff', '<f8'), ('lambda_r', '<f8'), ('vsig', '<f8'),
    ...  ('sfr01', '<f8'), ('sfr05', '<f8'), ('sfr1', '<f8'), ('sfr_area', '<f8'),
    ...  ('lgas', '<f8', (3,)), ('pos', '<f8', (3,)), ('vel', '<f8', (3,)),
    ...  ('lvec', '<f8', (3,)), ('nvec', '<f8', (3,))]

    lvec : spin direction
    nvec : normalized lvec


    >>> merger.dtype
    ... ['nstep', '<i8'), ('nout', '<i8'), ('id', '<i8'), ('idx', '<i8'),
    ... ('dist', '<f8'), ('orbitang', '<f8'), ('m', '<f8'),
    ... ('rel_pos', '<f8', (3,)), ('rel_vel', '<f8', (3,)),
    ... ('jorbit', '<f8', (3,)), ('spinang', '<f8'), ('spinmag', '<f8'),
    ... ('mstar', '<f8'), ('m_frac', '<f8'), ('m_gas', '<f8'), ('m_gas_cold', '<f8'),
    ... ('rgal', '<f8'), ('rgal_sum', '<f8'), ('reff_s', '<f8')]

    jorbit: cross(rel_pos, rel_vel)

    1. No need to interpolate :
        nstep, nout, id, idx, dist, rel_pos, rel_vel, jorbit

    2. Need to interpolate... well :
        orbitang, spinang, spinmag, m_frac,

    3. Simpler linear interpolation is fine :
        mstar, m_gas, m_gas_cold, size_s,   reff_s

     m_frac, size_p, rgal_sum, rgal_sum

    NOTE
    ----

    1. If sat_data["nvec"] has leading/trailing nan, the measurement at the point is ignored,
    but the nout point is included in the new_nout, at which more crude interpolation is done.
    Maybe that's OK as I am not interested in orbit-related values at the very first or last moment.

    2. Do I just throw away all the rest of the satellite properties??  - Maybe

    3. Merger data might be removed after cal_merger_epoch.


    4. Interpolate individual vector first, and then calculate angles
       so that I can normalize interpolated vectors beforehand.

    """

    tt = this_gal.finedata
    merger["lbt"] = nnza_all.a2b(merger["nout"], "nout", "lbt")
    i_nn = np.where(merger["nout"] >= this_gal.nout_data.min())[0]
    new_nout = merger["nout"][i_nn]
    i_mt = mtc.match_list_ind(tt["nout"], new_nout)

    new_lbt = merger["lbt"][i_nn]
    ###################################################
    # The rest: m_frac, rgal_sum
    # host mass, host size are required.
    merger["mstar"][i_nn] = nsm.interp_np(sat_data["lbt"], sat_data["mstar"],
                                      new_lbt, too_many_nan_frac=0.4)

    if do_smooth:
        #merger["mstar"] = nsm.smooth(merger["mstar"])
        merger["mstar"] = savgol_filter(merger["mstar"],
                                        win_size_large, 2, deriv=0, delta=1.0, axis=0, mode='interp', cval=0.0)
        merger["rel_pos"] = savgol_filter(merger["rel_pos"],
                                        win_size_large, 2, deriv=0, delta=1.0, axis=0, mode='interp', cval=0.0)
        merger["rel_vel"] = savgol_filter(merger["rel_vel"],
                                        win_size_small, 2, deriv=0, delta=1.0, axis=0, mode='interp', cval=0.0)

    merger["jorbit"]=np.cross(merger["rel_pos"], merger["rel_vel"])
    merger["dist"]=np.sqrt(np.einsum('...i,...i', merger["rel_pos"],merger["rel_pos"])) * 1e3 # in Kpc

    # There's no point interpolating points earlier than first main galaxy measurement
    # with no simple way to extrapolate different types of quantities reasonably.
    # Just don't do that.
    # Likewise, ignore nan / bad measurements.

    # But we have full access to main measurements.
    # Only new_nout is limited by sat length.
    #main_nvec = np.vstack(interp_vec_arr_to_arrays(this_gal.main_data["lbt"],
    #                                               this_gal.main_data["nvec"],
    #                                               new_lbt,
    #                                               normalize=True)).T
    main_nvec = tt["nvec"][i_mt]
    # j_orbit is not normalized when first added.
    # I will keep it unnormalized because the magnitude matters later.
    # So normalize it here, temporarily.

    # maintree & sattree
    merger["orbitang"][i_nn] = 180./np.pi*np.arccos(np.einsum('...i,...i',
                                                        main_nvec,
                                                        nsm.norm_3d(merger["jorbit"][i_nn])))
    try:
        sat_nvec_fine = np.vstack(nsm.interp_vec_arr_to_arrays(sat_data["lbt"],
                                                       sat_data["nvec"],
                                                       new_lbt,
                                                       normalize=True)).T
    except:
        return False

    #print(tt["mstar"][i_mt])
    #print(tt["lvec"][i_mt,0])
    #print(sat_nvec_fine[:,0])
    #print(main_nvec)
    merger["spinang"][i_nn] = 180./np.pi*np.arccos(np.einsum('...i,...i',
                                                   main_nvec,
                                                   sat_nvec_fine))

    ###################################################
    # Simpler(Linear) interpolations.
    # "mstar", "m_gas", "m_gas_cold", "size_s", "reff_s"
    for f_org, f_dest in zip(fields_to_interpol_from,fields_to_interpol_to):
        merger[f_dest][i_nn] = nsm.interp_np(sat_data["lbt"],
                                        sat_data[f_org],
                                        new_lbt,
                                        too_many_nan_frac=0.4)

    merger["m_frac"][i_nn] = merger["mstar"][i_nn]/nsm.smooth(this_gal.finedata["mstar"][i_mt])
    merger["rgal_sum"][i_nn] = nsm.smooth(merger["rgal"][i_nn]+this_gal.finedata["rgal"][i_mt])
    #
    # add gas frac? probably!
    #
    lvec_fine=np.vstack(nsm.interp_vec_arr_to_arrays(sat_data["lbt"],
                                                 sat_data["lvec"],
                                                 new_lbt,
                                                 normalize=False)).T
    merger["spinmag"][i_nn]=np.sqrt(np.einsum("...i,...i", lvec_fine, lvec_fine))

    return True # good = True


def add_results_to_gal(gals, nouts, prg_dir, out_base='./', verbose=False):
    out_dir = out_base + "lambda_results/"
    all_sample_ids= pickle.load(open(prg_dir + "all_sample_ids.pickle", "rb"))
    all_sample_idxs=pickle.load(open(prg_dir + "all_sample_idxs.pickle", "rb"))

    for nout in nouts:
        try:
            results_thisnout = pickle.load(open(out_dir+"results_{}.pickle".format(nout), "rb"))
        except:
            print("No pickle...")
            results_thisnout = []
            fn_all = glob(out_dir+"{}/result_sub_sample_{}_*.pickle".format(nout, nout))
            for fn in fn_all:
                idsnow = np.array(all_sample_ids[str(nout)])
                idxsnow = np.array(all_sample_idxs[str(nout)])
                # Some results have right idx, some are wrong...
                this_result = pickle.load(open(fn, "rb"))
                results_thisnout.extend(this_result)
            pickle.dump(results_thisnout, open(out_dir+"results_{}.pickle".format(nout), "wb"))
        # negative idx = phantom. So, what should I do with that?

        if nout == 782: print("There are {} measurements".format(len(results_thisnout)))

        good=0
        bad=0
        allids = np.array([agal.id for agal in results_thisnout])
        for this_gal in gals:
            # Main galaxy
            ind=np.where(this_gal.main_data["nout"]==nout)[0][0]
            itree = np.where(this_gal.finedata["nout"] == nout)[0][0]
            i_data = np.where(allids == this_gal.finedata[itree]["id"])[0]
            if len(i_data) == 1:
                data = results_thisnout[i_data[0]]
                nsm.add_main(this_gal.main_data, data, ind)

            for i, (merger, sat_data) in enumerate(zip(this_gal.mergers, this_gal.sat_data)):
                i_offset = np.where(merger["nout"] == nout)[0]
                if len(i_offset) > 0 :
                    try:
                        data = results_thisnout[np.where(allids == merger[i_offset]["id"])[0][0]]
                        good+=1
                    except:
                        bad+=1
                        continue
                        # What if there is no match?
                    ind=np.where(sat_data["nout"] ==nout)[0][0]
                    sat_data[ind]["mstar"] = data.mstar
                    sat_data[ind]["pos"] = (data.xc, data.yc, data.zc)
                    sat_data[ind]["vel"] = (data.vxc, data.vyc, data.vzc)

                    if hasattr(data, "rgal"):
                        sat_data[ind]["rgal"] = data.rgal
                        sat_data[ind]["reff"] = data.reff
                        sat_data[ind]["nvec"] = data.nvec
                        sat_data[ind]["lvec"] = data.lvec
                        if hasattr(data, "gas_results"):
                            sat_data[ind]["mgas"] = data.gas_results["mgas_tot"]
                            sat_data[ind]["mgas_cold"] = data.gas_results["mgas_cold"]
                        else:
                            sat_data[ind]["mgas"] = np.nan
                            sat_data[ind]["mgas_cold"] = np.nan
                    else:
                        sat_data[ind]["rgal"] = np.nan
                        sat_data[ind]["reff"] = np.nan
                        sat_data[ind]["nvec"] = np.nan
                        sat_data[ind]["lvec"] = np.nan
                        sat_data[ind]["mgas"] = np.nan
                        sat_data[ind]["mgas_cold"] = np.nan

        if verbose: print(good, bad)





def find_merger_ini_ends(j_orbit_mag, k, j0, lbt, lbt_min_j_decrease,
                         j0_percentile=80,
                         dt_merger_ini = 0.5,
                         dt_merger_fi=-0.2,
                         ii_dt_min = 20):
    """
        Parameters
        ----------
        dt_merger_fi:
            what is this...?

        NOTE
        ----
        Next merger can happen 0.75 * lbt_min_j_decrease after the end of earlier merger.
        If not, two event are merged to be one longer merger.



    """
    i_merger_ini=[]
    i_merger_end=[]
    n_min_decrease_j = np.argmax(lbt > lbt[k] - lbt_min_j_decrease)

    seek_ini=True
    while k > 0:
        # j should decrease for at least n_min_decrease (0.5Gyr).
        # As long as j does not grow above j0.

        if seek_ini:
            # np.all(j_orbit_mag[k-n_min_decrease_j:k] < j0)
            # This condition requires j_orbit_mag can not fluctuate to get shot up above j0.
            # Isn't it too strict?
            if j_orbit_mag[k] >= j0 and np.all(j_orbit_mag[k-n_min_decrease_j:k] < j0):

                #print("New j0 {:.2f} at ind = {}".format(j0,k))
                if len(i_merger_end) > 0:
                    # If the next merger starts before the end of the previous merger
                    # or if the next one starts very soon, then merge them.
                    #if i_merger_end[-1] <= i+min_dstep_mergers:
                    #
                    if i_merger_end[-1] <= k+n_min_decrease_j * 0.75:
                        #print("remove", i_merger_end[-1])
                        i_merger_end.pop(-1)
                        seek_ini=False
                        continue
                        #print("not continued")
                i_merger_ini.append(k)
                k -= n_min_decrease_j # Because this far has been tested.
                seek_ini=False
                #print("modified k",k)
        else:
            # if it rise again, above the initial j0
            if j_orbit_mag[k] <= j0 and j_orbit_mag[k-1] > j0:
                #print("end", k)
                i_merger_end.append(np.argmax(lbt - lbt[k] > dt_merger_fi))
                seek_ini=True

                if len(i_merger_ini) < len(i_merger_ini_all):
                    i_ini = np.argmax(i_merger_ini_all[:-len(i_merger_end)] < i_merger_end[-1])
                    #print("new i_ini", i_ini)
                    ii_dt = np.argmax(lbt - lbt[i_ini] > dt_merger_ini)
                    new_ini = min([max([ii_dt, i_ini + ii_dt_min]), len(j_orbit_mag) -1])
                    # minimum j0 = j0 at 3R
                    #j0 = np.max([np.median(j_orbit_mag[i_ini:new_ini]) * j0_merger_frac, j_orbit_mag[i_ini]])
                    j0 = np.percentile(j_orbit_mag[i_ini:new_ini], j0_percentile)
                    k =i_ini +1
                    # update n_min_decrease_j
                    n_min_decrease_j = np.argmax(lbt > lbt[k] - lbt_min_j_decrease)
        k -=1

    if not seek_ini:
        i_merger_end.append(0)

    return np.array(i_merger_ini), np.array(i_merger_end)
