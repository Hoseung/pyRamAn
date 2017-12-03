from utils import hagn
import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt 
import galaxymodule  # needed for result_sub_sample_**.pickle
import utils.match as mtc 
import pickle
from utils import hagn
import os
import tree.halomodule as hmo 
import cell_chunk_module as ccm 
from utils import cosmology
from load.info import Info
#import P_measure_fine
from rot2 import density_measure as denm
from time import time
import importlib
from glob import glob
from rot2 import new_serial as ns
import rot2.new_serial_modules as nsm
from scipy.signal import savgol_filter


def gather_sample_gals(all_final_idxs, nnza_cell, nnza_all, prg_dir,istep_max,
                        nstep_main_too_short = 500):

    gals=[]
    for j, this_idx in enumerate(all_final_idxs):
        fname = prg_dir + "{}_adp.pickle".format(this_idx)
        if not os.path.isfile(fname):
            # dump_prgs broken
            print(j,"Missing prgs file for", this_idx)
            continue
        adp = pickle.load(open(fname, "rb"))

        # abort if the tree is too short (in all_cell nstep).
        if min(adp[0][0]["nstep"][adp[0][0]["nstep"] > 0]) > nstep_main_too_short:
            print("skip")
            continue
        # Append age to maintree and mainresult.
        this_gal = ns.Maingal(adp, nnza_all, istep_max, nnza_cell)
        this_gal.add_mergers_tree
        #cnt_merger=0
        bad_main=False
        gals.append(this_gal)

    return gals


def nout_results_to_gal_serial(gals, nouts, out_base):
    out_dir = out_base + "lambda_results/"

    for nout in nouts:
        # All results at this nout in a list.
        print(nout)
        try:
            results_thisnout = pickle.load(open(out_dir+"results_{}.pickle".format(nout), "rb"))
        except:
        # Get right IDx
            print("No pickle...")
            results_thisnout = []
            fn_all = glob(out_dir+"{}/result_sub_sample_{}_*.pickle".format(nout, nout))
            for fn in fn_all:
                # Some results have right idx, some are wrong...
                this_result = pickle.load(open(fn, "rb"))
                results_thisnout.extend(this_result)
            pickle.dump(results_thisnout, open(out_dir+"results_{}.pickle".format(nout), "wb"))
        # negative idx = phantom. So, what should I do with that?

        if nout == 782: print("There are {} measurements".format(len(results_thisnout)))

        # Distribute results to each gal.
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
                    # put values.
                    #main_tree = this_gal.finedata[mtc.match_list_ind(this_gal.finedata["nout"], merger["nout"])]
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
                            sat_data[ind]["mgas"] = 0
                            sat_data[ind]["mgas_cold"] = 0
                    else:
                        sat_data[ind]["rgal"] = np.nan
                        sat_data[ind]["reff"] = np.nan
                        sat_data[ind]["nvec"] = np.nan
                        sat_data[ind]["lvec"] = np.nan
                        sat_data[ind]["mgas"] = np.nan
                        sat_data[ind]["mgas_cold"] = np.nan

        print(good, bad)


def sat_data_to_fine(this_gal, merger, sat_data, nnza_all,
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
    fields_to_interpol_from=["mstar", "rgal",  "reff"] # ignore gas for the moment
    fields_to_interpol_to  =["mstar", "rgal", "reff_s"]

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


def find_merger_ini_ends(j_orbit_mag, k, j0, lbt, lbt_min_j_decrease, i_merger_ini_all,
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
        if i_merger_ini[-1] < 5:
            i_merger_ini.pop(-1)
        if len(i_merger_ini) > len(i_merger_end):
            i_merger_end.append(0)            
    
    return np.array(i_merger_ini), np.array(i_merger_end)



def cal_merger_props(gals, nouts_all,
                     verbose=False,
                     j_smooth_w = 51,
                     smooth_poly_deg = 5,
                     rscale = 3.0,
                     j0_mean_window = 0.5, # take median of j over this time span
                     # Merger criterion by angular momentum.
                     # Emperically determined..
                     j0_merger_frac = 1.0,
                     j0_percentile = 0.8, # top 20% of j within the window.
                     lbt_min_j_decrease = 0.5):

    """
    parameters
    ----------
    j_smooth_w = 51 
        Very large smoothing window. 

    Note one removing lis element.

    I want to remove mergers that do not reach close enough to the host.
    Both this_gal.mergers[i] and this_gal.sat_data[i] need to be removed 
    while the loop runs through this_gal.mergers. 


    Merging sat always within the distance criterion
    - If it's a later fraction of merging sat tree, 
    a merging sat could stay within the distance criterion for all the time. 
    This slips away my criterion that a sat should come from outer than a certain distance. 
    Should I include them as a healthy merger? or just ignore them?
    Probably the earlier part of the tree is in the this_gal.mergers list ... 
    => Remove them. Or, I can try fixing it in earlier steps. 


    """

    dl_all = []
    dj_all = []

    j_smooth_w = 51 
    smooth_poly_deg = 5
    rscale = 3.0
    j0_mean_window = 0.5 # take median of j over this time span
    # Merger criterion by angular momentum. 
    # Emperically determined.. 
    j0_merger_frac = 1.0
    j0_percentile = 0.8 # top 20% of j within the window.
    #

    lbt_min_j_decrease = 0.5

    for igal, tg in enumerate(gals):
        tg.merger_meta = {"i_ini":[],
                          "i_end":[],
                          "m_ratio":[],
                          "dj":[],
                          "dl":[],
                          "m_ratio_ratio":[],
                          "frac_merger_time":0}
        
        bad_mergers=[]    
        
        for imerger, merger in enumerate(tg.mergers):
            # global meta
            #ind_offset = np.where(tg.finedata["nout"] == merger["nout"][0])[0][0]
            ind_offset = tg.finedata["nout"][0] - merger["nout"][0]

            # global data
            lbt = merger["lbt"]
            m_ratio = merger["m_frac"]

            r_dist_frac = merger["dist"]/merger["rgal_sum"]

            # find j_orbit_max
            #print(len(merger))
            if len(merger) == 15:
                j_smooth_w = 13
            else:
                j_smooth_w = 15
            smooth_jorb = savgol_filter(merger["jorbit"],
                                        j_smooth_w,
                                        smooth_poly_deg,
                                        deriv=0, delta=1.0, axis=0, mode='interp', cval=0.0)
            j_orbit_mag = np.sqrt(np.einsum("...i,...i", smooth_jorb, smooth_jorb))

            # merger_init candidates
            # When dist crosses the threshold
            i_merger_ini_all = np.where((r_dist_frac[1:] > rscale) * (r_dist_frac[:-1] < rscale))[0]
            if len(i_merger_ini_all) > 0:
                #print("All ini", i_merger_ini_all)
                i_ini_ini = i_merger_ini_all[-1] # Earliest merger ini

                # If there is not enough tree to take mean from, 
                # just take as much as possible.
                if (lbt[-1] - lbt[i_ini_ini]) < j0_mean_window:
                    if (lbt[-1] - lbt[i_ini_ini]) > 0:
                        ii_dt=len(lbt)                    
                else:
                    ii_dt = np.argmax(lbt - lbt[i_ini_ini] > j0_mean_window)
                
                new_ini = min([ii_dt, len(j_orbit_mag) -1])

                # median jorbit at the initial time over a certain window. 
                j0 = np.percentile(j_orbit_mag[i_ini_ini:new_ini], j0_percentile)

                # Look for merger beginning and end.
                i_merger_ini, i_merger_end = find_merger_ini_ends(j_orbit_mag, new_ini, j0, lbt, lbt_min_j_decrease, i_merger_ini_all)
                
            else:
                if verbose: print("No mergers by_distance... SKIP")
                bad_mergers.append(imerger)
                #nsm.check_merger_prps(merger, tg.fidx, merger["idx"][0])
                continue

            if len(i_merger_ini) == 0:
                if verbose: print("No mergers by_j... SKIP")
                bad_mergers.append(imerger)
                # remove the merger AFTER the loop is over.            
                continue
                
            # Merger meta data 
            tg.merger_meta["i_ini"].append(i_merger_ini + ind_offset)
            tg.merger_meta["i_end"].append(max([0,i_merger_end + ind_offset])) # ends at 0.
            tg.merger_meta["m_ratio"].append(m_ratio[i_merger_ini])
            tg.merger_meta["dj"].append(j_orbit_mag[i_merger_ini] - j_orbit_mag[i_merger_end])
            
        #plt.savefig("merger_pos_test_{}_{}.png".format(igal, tg.fidx))
        
        # Loop is done. Remove unnecessary mergers.
        # Should remove from larger index.
        for i_bad in np.sort(bad_mergers)[::-1]:
            tg.mergers.pop(i_bad)
            tg.sat_data.pop(i_bad)
        
        # weights of each merger by merger mass ratio.

        #main_d_lambda = tg.finedata["lambda_r"][:-1] - tg.finedata["lambda_r"][1:]
        all_weights = np.zeros(len(nouts_all) + 20) # Why 20?
        #n_major_arr = np.zeros(len(tg.finedata))
        #n_minor_arr = np.zeros(len(tg.finedata))

        #tg.n_major=0
        #tg.n_minor=0

        #for merger in tg.mergers:
        #    #if hasattr(merger, "i_ini"):
        #    for ini, end, mratio in zip(merger.i_ini, merger.i_end, merger.m_ratio):
        #        #print(ini, end)
        #        all_weights[end:ini+1] += mratio
        #        if mratio > mratio_Major:
        #            n_major_arr[end:ini+1] += 1
        #            tg.n_major+=1
        #        elif mratio > mratio_minor:
        #            n_minor_arr[end:ini+1] += 1
        #            tg.n_minor+=1

        # Merger history types.
        #tg.n_major_arr=n_major_arr
        #tg.n_minor_arr=n_minor_arr


        # A window to take the mean of delta lambda?
        # Not needed, I've already smoothed the values.     
        mm = tg.merger_meta
        if len(mm["i_ini"]) > 0:
            # Compute allweights first
            for ini, end, mr in zip(mm["i_ini"], mm["i_end"], mm["m_ratio"]):
                for i_ini, i_end, mratio in zip(ini, end, mr):
                    all_weights[i_end:i_ini] += mratio

            for ini, end, mr, tdj in zip(mm["i_ini"], mm["i_end"], mm["m_ratio"], mm["dj"]):
                for i_ini, i_end, mratio, thisdj in zip(ini, end, mr, tdj):
                    thisdl = tg.finedata["lambda_r"][i_end:i_ini] - tg.finedata["lambda_r"][i_end+1:i_ini+1]
                    thismr = mratio/all_weights[i_end:i_ini]
                    mm["m_ratio_ratio"].append(thismr)
                    mm["dl"].append(np.sum(thisdl * thismr))
                    #dl_all.append(np.sum(thisdl * thismr))
                    #dj_all.append(thisdj)
        mm["frac_merger_time"] = np.sum(all_weights > 0) / len(all_weights)    
                    
        # return dl_all, dj_all
        #plt.close("all")



def measure_density(gals, nnza_all, nnza_cell,
                    nouts, nouts_all, sim_base):
    for this_gal in gals:
        this_gal.env=np.zeros(len(nnza_all.nnza), dtype=[("d5", "<f8"),
                                                         ("d10", "<f8"),
                                                         ("d50", "<f8"),
                                                         ("P_tidal", "<f8")])

        this_gal.env_short=np.zeros(len(nnza_cell.nnza), 
                                                  dtype=[("P_tidal_h", "<f8"),
                                                         ("host_d_id", "<i4"),
                                                         ("host_d_m", "<f8"),
                                                         ("host_t1r_id", "<i4"),
                                                         ("host_t1r_m", "<f8"),
                                                         ("host_t2r_id", "<i4"),
                                                         ("host_t2r_m", "<f8"),
                                                         ("pos1", "<f8", (3,)),
                                                         ("Rvir1", "<f8"),
                                                         ("r_dist1", "<f8"),
                                                         ("pos2", "<f8", (3,)),
                                                         ("Rvir2", "<f8"),
                                                         ("r_dist2", "<f8")])

    """ 
        this_gal.main_data will be removed.
        Use finedata whenever possible.
    """
    for i,nout in enumerate(nouts_all):
        if nout in [584,585,359,293,294]:
            continue

        #if nout not in nouts:
        #    continue
        gdata = pickle.load(open(sim_base+"GalaxyMaker/gal_pickle/gcat_{}.pickle".format(nout),"rb"))
        print("Now ", nout)
        info = Info(base=sim_base, nout=nout)
        denm.density_D2N(gdata, info, gals, Ns=[5, 10, 50])
        dt_fine = nnza_all.nnza["lbt"][i]-nnza_all.nnza["lbt"][i-1]
        denm.measure_P(gdata, info, gals, dt_fine, short=False)


        # Only 63 snapshots
        if nout not in nnza_cell.nnza["nout"]:
            continue
        else:
            hdata = pickle.load(open(sim_base+"halo/DM_pickle/hcat_{}.pickle".format(nout),"rb"))
            inout_cell = np.where(nnza_cell.nnza["nout"] == nout)[0]
            dt = nnza_cell.nnza["lbt"][inout_cell-1]-nnza_all.nnza["lbt"][inout_cell]
            denm.measure_P(hdata, info, gals, dt, short=True)

            all_ids_now=[]
            sample_inds=[]
            for ii, tg in enumerate(gals):
                igal_now = np.where(tg.finedata["nout"]==nout)[0]
                #print(igal_now)
                if len(igal_now) > 0:
                    all_ids_now.extend(tg.finedata["id"][igal_now])
                    sample_inds.append(ii)
            sample_inds = np.array(sample_inds)
            all_ids_now = np.array(all_ids_now)

            i_cell = np.where(nouts == nout)[0]
            direct_hosts, largest_hosts, largest_hosts2 = denm.match_halo_gal(all_ids_now, gdata, hdata)
            for j, igal in enumerate(sample_inds):
                this_gal = gals[igal]
                this_gal.env_short["host_d_id"][i_cell]   = direct_hosts[j]["id"]
                this_gal.env_short["host_d_m"][i_cell]    = direct_hosts[j]["mvir"]
                this_gal.env_short["host_t1r_id"][i_cell] = largest_hosts[j]["id"]
                this_gal.env_short["host_t1r_m"][i_cell]  = largest_hosts[j]["mvir"]
                this_gal.env_short["host_t2r_id"][i_cell] = largest_hosts2[j]["id"]
                this_gal.env_short["host_t2r_m"][i_cell]  = largest_hosts2[j]["mvir"]

                this_gal.env_short["pos1"][i_cell] = (largest_hosts[j]["x"]-0.5,
                                                      largest_hosts[j]["y"]-0.5,
                                                      largest_hosts[j]["z"]-0.5)
                this_gal.env_short["pos1"][i_cell]*=  info.pboxsize
                this_gal.env_short["Rvir1"][i_cell] =  largest_hosts[j]["rvir"] * info.pboxsize

                this_gal.env_short["pos2"][i_cell] = (largest_hosts2[j]["x"]-0.5,
                                                      largest_hosts2[j]["y"]-0.5,
                                                      largest_hosts2[j]["z"]-0.5)
                this_gal.env_short["pos2"][i_cell]*=info.pboxsize
                this_gal.env_short["Rvir2"][i_cell] = largest_hosts2[j]["rvir"] * info.pboxsize


