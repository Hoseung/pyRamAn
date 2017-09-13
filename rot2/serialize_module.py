import numpy.lib.recfunctions as recf
import pickle
import os
import numpy as np
from rot2 import serialize_results
from glob import glob

def idxfromid(org_ids, ids_ref, idxs_ref):
    i_sort = np.argsort(ids_ref)
    ids_ref_new = ids_ref[i_sort]
    idxs_ref_new = idxs_ref[i_sort]

    return idxs_ref_new[np.argsort(org_ids).argsort()]


def get_all_results(nouts,
                    prg_dir = "./test_direct_prgs_gal/",
                    out_dir = "./lambda_results/"):
    # ALL ALL results.
    all_sample_ids=pickle.load(open(prg_dir + "all_sample_ids.pickle", "rb"))
    all_sample_idxs=pickle.load(open(prg_dir + "all_sample_idxs.pickle", "rb"))
    allresults=[]
    for nout in nouts:
        allresults_thisnout = []
        # Better that the list of sample is stored in a separate file, instead of
        # trying to read ALL files in a directory...
        fn_all = glob(out_dir+"{}/result_sub_sample_{}_*.pickle".format(nout, nout))
        # Get right IDx
        for fn in fn_all:
            idsnow = np.array(all_sample_ids[str(nout)])
            idxsnow = np.array(all_sample_idxs[str(nout)])

            # Some results have right idx, some are wrong...
            this_result = pickle.load(open(fn, "rb"))
            allidxs = np.array([agal.idx for agal in this_result])
            if max(allidxs) < 1e6:
                allidxs = idxsnow[mtc.match_list_ind(idsnow, allidxs)]
                #idxfromid(allidxs, idsnow, idxsnow)
                for idx, agal in zip(allidxs, this_result):
                    agal.idx =idx
            allresults_thisnout.extend(this_result)

        allresults.append(allresults_thisnout)
    return allresults


def smooth(x, beta=5, window_len=20, monotonic=False, clip_tail_zeros=True):
    """
    kaiser window smoothing.

    If len(x) < window_len, window_len is overwritten to be len(x).
    This ensures to return valid length fo array, but with modified window size.

    Parameters
    ----------
        window_len = 20

        monotoinc =

        clip_tail_zereos = True
            returned array is shorter than the original.

    beta = 5 : Similar to Hamming


    """
    lx = len(x)
    if clip_tail_zeros:
        x = x[:max(np.where(x > 0)[0])+1]

    if monotonic:
        """
        if there is an overall slope, smoothing may result in offset.
        compensate for that.
        """
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y=np.arange(len(x)))
        xx = np.arange(len(x)) * slope + intercept
        x = x - xx

    # extending the data at beginning and at the end
    # to apply the window at the borders
    window_len = min([window_len, len(x)])
    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]] # concatenate along 0-th axis.
    # periodic boundary.
    w = np.kaiser(window_len,beta)
    y = np.convolve(w/w.sum(), s, mode='valid')
    aa = len(y)-lx
    if monotonic:
        return y[int(window_len/2):len(y)-int(window_len/2) + 1] + xx
    else:
        #return y[int(window_len-1/2)-2:len(y)-int((window_len-1)/2)]
        return y[int(aa/2):int(aa/2)+len(x)]


def cal_merger_epoch(serial_results, tc,
                 nnza_cell, nnza_all,
                 prg_dir, serial_out_dir,
                 nstep_too_short_main = 100,
                 do_plot=False,
                 merger_something = False,
                 remove_small = True,
                 merger_criterion="by_j",
                 rscale=3.0,
                 dt_merger_ini = 0.5,
                 dt_merger_fi  = -0.2,
                 min_dstep_mergers = 15,
                 ii_dt_min = 20,
                 n_min_decrease_j = 15,
                 m_ratio_min = 1/50,
                 j0_merger_frac=1.0,
                 j_smooth_window = 25,
                 mratio_Major=0.25,
                 mratio_minor=0.02):
    """
        Parameters
        ----------
        min_dstep_mergers : 15
        if a merger starts within min_dstep_mergers after the end of the last merger, merge them.
        ii_dt_min : 20
        j0 is the median of j_orbit between the beginning of the merger and at least ii_dt_min nsteps earlier.
        n_min_decrease_j : 15
        for a merger to start, j0 must decrease for at least for this long.
    """

    # Determine beginning and ending of mergers.
    kpc_in_km = 3.08567758e16
    yr_in_s = 3600*24*365
    #scale_cross = 3.0

    #merger_criterion = ["by_dist", "by_j"][1]
    d_frac_min=[]
    delta_Msat=[]

    # j, orbital angular momentum is the central parameter in determining the merger epoch (beginning and end).
    # For a more reliable determination, smooth it quite a bit.
    if do_plot:
        fig, axs = plt.subplots(2,2, sharex=True)
        fig.set_size_inches(8,6)
        axs = axs.ravel()

    dl_all = []
    dj_all = []

    Major_dominated=[] # Major merger overlaps totally
    minor_dominated=[] #
    mixed_mergers=[]
    rare_mergers=[]


    for i, this_gal in enumerate(serial_results[7:19]):
        try:
            len(this_gal.data[0][0]) > 1
        except:
            continue
        print(i, "IDX=",this_gal.fidx, end="\r")

        lbt_main_gal = nnza_all.a2b(this_gal.finearr["nstep"], "nstep", "lbt")

        this_gal.n_major=0
        this_gal.n_minor=0

        this_gal.n_major_arr=np.zeros(len(this_gal.finearr))
        this_gal.n_minor_arr=np.zeros(len(this_gal.finearr))


        for merger in this_gal.mergers[21:]:
            zred_fine = nnza_all.a2b(merger.main_tree["nstep"],"nstep","zred")
            zred_sat = nnza_cell.a2b(merger.sat["nstep"],"nstep","zred")
            i_lbt_end = np.argmax(zred_fine > zred_sat.max())

            min_nstep_sat = merger.main_tree["nstep"]
            merger.sattree = merger.sattree[i_lbt_end]
            merger.main_tree = merger.main_tree[i_lbt_end]
            main_tree_part = merger.main_tree
            # dist
            rel_pos = main_tree_part["pos"] - merger.sattree["xp"]
            dist3d  = np.sqrt(np.einsum("...i,...i",rel_pos,rel_pos)) * 1e3 # in kpc
            if np.sum(dist3d > 0) < 5:
                continue

            main_nstep_ini = main_tree_part["nstep"][-1]
            main_nstep_fi  = main_tree_part["nstep"][0]
            main_finearr = this_gal.finearr[(this_gal.finearr["nstep"] <=main_nstep_fi) *\
                                              (this_gal.finearr["nstep"] >=main_nstep_ini)]

            #print("mtp",main_tree_part["nstep"], "mf",main_finearr["nstep"])
            #print("mst",merger.sattree["nstep"], merger.sat["nstep"])
            #rel_pos = main_finearr["pos"] - merger.sat
            # vel
            rel_vel=main_tree_part["vel"] - merger.sattree["vp"]
            v3d = np.sqrt(np.einsum("...i,...i", rel_vel, rel_vel))

            m_ratio = merger.sattree["m"]/main_tree_part["mstar"]
            if max(m_ratio) < m_ratio_min:
                print("Too small", max(m_ratio[-50:]))
                if remove_small:
                    this_gal.mergers.remove(merger)
                continue
            # galaxy radius over fine nstep.
            lbt = tc.zred2gyr(nnza_all.a2b(main_tree_part["nstep"],"nstep","zred"), z_now=0)


            rp_fine = main_finearr["rgal"]
            r_s_smoothed = smooth(merger.sat["rgal"], window_len=5, clip_tail_zeros=False)
            fine_nout = main_tree_part["nstep"]+30 # max(nout) - max(nstep)
            rs_fine = np.interp(fine_nout[::-1], merger.sat["nout"][::-1], r_s_smoothed[::-1])[::-1]

            # Gas mass
            #mg_s = smooth(merger.sat["mgas"], window_len=5, clip_tail_zeros=False)
            #mgc_s = smooth(merger.sat["mgas_cold"], window_len=5, clip_tail_zeros=False)
            mg_s_fine = np.interp(fine_nout[::-1], merger.sat["nout"][::-1], merger.sat["mgas"][::-1])[::-1]
            mgc_s_fine = np.interp(fine_nout[::-1], merger.sat["nout"][::-1], merger.sat["mgas_cold"][::-1])[::-1]

            # Arrays : mass, sfr, vsig, j_orbit
            # scalars : Mgas,
            s_mass = np.interp(fine_nout[::-1], merger.sat["nout"][::-1], merger.sat["mstar"][::-1])[::-1]
            s_sfr  = np.interp(fine_nout[::-1], merger.sat["nout"][::-1], merger.sat["sfr01"][::-1])[::-1]
            s_vsig = np.interp(fine_nout[::-1], merger.sat["nout"][::-1], merger.sat["vsig"][::-1])[::-1]
            #print(merger.sat["mstar"])
            #print(s_mass, main_finearr["mstar"])
            m_ratio = s_mass / main_finearr["mstar"]

            r_tot = rp_fine + rs_fine

            # cumulative min dist Vs dM
            d_frac_min.extend(np.minimum.accumulate((dist3d/r_tot)[-2::-1]))
            delta_Msat.extend((s_mass[:-1]-s_mass[1:])/s_mass[1:])


            #v_ini = v3d[i_ini]
            #print("Velocity at the beginning of the merger", v_ini)
            #print("Note that the relative velocity depends A LOT on when it is measured/observed.")
            # kms -> kpc_per_Gyr
            #t_cross = scale_cross*rgal_ini*kpc_in_km/(v_ini*31536000*1e9)
            # End of merger
            #print("Rgal = {}, Crossing time {} Gyr".format(rgal_ini, t_cross))
            #i_fin = np.argmax(lbt > (t_ini + t_cross)) #  ???

            j_orbit = np.cross(rel_pos,rel_vel)
            #print(m_ratio, j_orbit, ":::::")
            j_orbit_mag = np.sqrt(np.einsum("...i,...i", j_orbit,j_orbit))#
            # Smoothing
            #print(j_orbit_mag.shape)
            j_orbit_mag = smooth(j_orbit_mag, window_len=min([j_smooth_window, int(sum(j_orbit_mag > 0)/2)]))
            #j_orbit_mag = j_orbit_mag[:len(m_ratio)]
            #print(j_orbit_mag.shape)
            if len(j_orbit_mag) < 5:
                print("Bad sat measurements, skip...")
                continue

            if len(merger.sat) == 1:
                spin_ang = np.repeat(merger.merger_arr["spinang"], len(merger.sattree))
                spin_mag = np.repeat(merger.merger_arr["spinmag"], len(merger.sattree))
                orbit_ang = np.repeat(merger.merger_arr["orbitang"], len(merger.sattree))
            else:# len(merger.sat) < 5:
                #smooth(merger.main_part["rgal"], window_len=5, clip_tail_zeros=False)
                spin_ang = np.interp(fine_nout[::-1], merger.sat["nout"][::-1], merger.merger_arr["spinang"])
                spin_mag = np.interp(fine_nout[::-1], merger.sat["nout"][::-1], merger.merger_arr["spinmag"])
                orbit_ang = np.interp(fine_nout[::-1], merger.sat["nout"][::-1], merger.merger_arr["orbitang"])

            # [:-1] Later
            # [1:] earlier
            #i_pericenter = # smoothed dist3d.
            i_j_max = 4 + np.argmax(j_orbit_mag[4:])
            # because j should decrease during the merger,
            # it should appear after the j_max.
            i_merger_ini_all = np.where((dist3d[1:i_j_max]  > rscale*(r_tot[1:i_j_max])) \
                               *(dist3d[:i_j_max-1] < rscale*(r_tot[:i_j_max-1])))[0]
            if len(i_merger_ini_all) > 0 and i_merger_ini_all[0] > 0:
                i_ini = i_merger_ini_all[-1]
                ii_dt = np.argmax(lbt - lbt[i_ini] > dt_merger_ini)
                #print(ii_dt, i_ini)
                new_ini = min([max([ii_dt, i_ini + ii_dt_min]), len(j_orbit_mag) -1])
                # minimum j0 = j0 at 3R
                j0 = np.max([np.median(j_orbit_mag[i_ini:new_ini]) * j0_merger_frac, j_orbit_mag[i_ini]])

                i_merger_ini=[]
                i_merger_end=[]
                seek_ini=True
                k = new_ini
                #print("new i_ini", i)
                #print(j_orbit_mag)
                while k > 0:#i in range(new_ini-1,0,-1):
                    if seek_ini:
                        if j_orbit_mag[k] >= j0 and np.all(j_orbit_mag[k-n_min_decrease_j:k] < j0):
                            #print(i, j0, n_min_decrease_j, j_orbit_mag[i-n_min_decrease_j:i])
                            #print("New j0", j0)
                            if len(i_merger_end) > 0:
                                #print(i, i_merger_end[-1])
                                # If the next merger starts before the end of the previous merger
                                # or if the next one starts very soon, then merge them.
                                if i_merger_end[-1] <= i+min_dstep_mergers:
                                    #print("remove", i_merger_end[-1])
                                    i_merger_end.pop(-1)
                                    seek_ini=False
                                    continue
                                    #print("not continued")
                            i_merger_ini.append(k)
                            k -= n_min_decrease_j # Because this far has been tested.
                            seek_ini=False
                    else:
                        if j_orbit_mag[k] <= j0 and j_orbit_mag[k-1] > j0:
                            #print("end", i)
                            i_merger_end.append(np.argmax(lbt - lbt[k] > dt_merger_fi))
                            seek_ini=True

                            if len(i_merger_ini) < len(i_merger_ini_all):
                                i_ini = np.argmax(i_merger_ini_all[:-len(i_merger_end)] < i_merger_end[-1])
                                #print("new i_ini", i_ini)
                                ii_dt = np.argmax(lbt - lbt[i_ini] > dt_merger_ini)
                                new_ini = min([max([ii_dt, i_ini + ii_dt_min]), len(j_orbit_mag) -1])
                                # minimum j0 = j0 at 3R
                                j0 = np.max([np.median(j_orbit_mag[i_ini:new_ini]) * j0_merger_frac, j_orbit_mag[i_ini]])
                                k =i_ini +1
                    k -=1

                if len(i_merger_ini) == 0:
                    print("BAD by_j... SKIP")
                    continue

                d_i_main_sat = this_gal.finearr["nstep"][0]-main_finearr["nstep"][0]
                # if ended without finding end of merger.
                # This is the eventual end.
                if not seek_ini:
                    #i_merger_end_main_gal =
                    # If there is no time where  lbt_main_gal - lbt[i] > dt_merger_fi
                    # because lbt[i] is close to z=0, np.argmax returns 0, which is perfectly fine.
                    i_merger_end.append(np.argmax(lbt_main_gal - lbt[0] > dt_merger_fi) - d_i_main_sat)
                    #i_merger_end.append(0)
                #print("Last merger", i_merger_end[-1])

                i_merger_ini = np.array(i_merger_ini)
                i_merger_end = np.array(i_merger_end)
                merger.i_ini = i_merger_ini + d_i_main_sat
                merger.i_end = i_merger_end + d_i_main_sat # to be compatible with this_gal.finetree
                # revert for pretty plot.
                i_merger_end[-1]=0
                #print(i, this_gal.fidx, i_merger_ini, i_merger_end, d_i_main_sat)
                merger.m_ratio = m_ratio[i_merger_ini] # sat property.
                merger.m_gas = mgc_s_fine[i_merger_ini] # sat property.
                #print(m_ratio.shape, j_orbit_mag.shape, [i_merger_end])
                #print((j_orbit_mag*m_ratio)[i_merger_ini])
                merger.dj = (j_orbit_mag*m_ratio)[i_merger_end]-(j_orbit_mag*m_ratio)[i_merger_ini]

                # IF there is a fly-by, separate it out.
                if do_plot:
                    j_orbit_mag *= m_ratio
                    # only for visualization purpose.
                    # In determining merger duration, no mass ratio is considered.
                    # Merger should be determined in its own right, without considering the effect on the primary galaxy.
                    #
                    #

                    axs[0].plot(merger.sattree["nstep"], dist3d)
                    #print(i_merger_ini, dist3d[i_merger_ini])
                    axs[0].scatter(merger.sattree["nstep"][i_merger_ini], dist3d[i_merger_ini], c="g",marker="+", s=50)
                    axs[0].scatter(merger.sattree["nstep"][i_merger_end], dist3d[i_merger_end], c="r", marker="^", s=50)
                    axs[0].plot(merger.sattree["nstep"], (r_tot)*rscale, "--")
                    axs[0].set_ylabel("Distance in kpc")
                    axs[0].set_xlabel("snapshot, smaller the later")
                    axs[0].set_ylim([0,1500])

                    axs[1].plot(this_gal.finearr["nstep"], this_gal.finearr["lambda_r"], color="g")
                    axs1 = axs[1].twinx()
                    axs1.plot(this_gal.finearr["nstep"], np.log10(this_gal.finearr["mgas"]+1)+1, color="r")
                    axs1.plot(this_gal.finearr["nstep"], np.log10(this_gal.finearr["mgas_cold"]+1)+1, color="b")
                    axs1.set_ylim([6,11])
                    axs[1].set_ylabel(r"$\lambda$")
                    axs[1].set_ylim([0,0.8])
                    axs[2].scatter(merger.sattree["nstep"], j_orbit_mag, c=orbit_ang, cmap="coolwarm", vmin=0, vmax=180)
                    axs[2].scatter(merger.sattree["nstep"][i_merger_ini], j_orbit_mag[i_merger_ini], c="g",marker="+", s=50)
                    axs[2].scatter(merger.sattree["nstep"][i_merger_end], j_orbit_mag[i_merger_end], c="r", marker="^", s=50)
                    axs[2].set_ylabel("j_obit")

                    #axs[3].plot(this_gal.finearr["nstep"], np.log10(this_gal.finearr["mstar"]+1))
                    # SFR
                    #young_stars=np.array([tg.sfr_results["hist"][0]*tg.sfr_results["area"] for tg in this_gal.data[0][0]])
                    #sfr_dt_gyr =np.array([tg.sfr_results["hist_dt"] for tg in this_gal.data[0][0]])
                    sfr_fine = this_gal.finearr["sfr01"]#np.interp(main_tree_part["time"], lbt_main, young_stars/sfr_dt_gyr)#*1e-9 # Msun in yr
                    #ssfr = young_stars*sfr_dt_gyr*1e9/this_gal.main_arr.mstar
                    ssfr = sfr_fine / this_gal.finearr["mstar"]
                    axs[3].scatter(this_gal.finearr["nstep"], np.log10(this_gal.finearr["mstar"]+1), c=ssfr)
                    axs[3].plot(merger.sattree["nstep"], np.log10(s_mass+1))
                    axs[3].plot(merger.sattree["nstep"], np.log10(mg_s_fine+1)+1, "--")
                    axs[3].plot(merger.sattree["nstep"], np.log10(mgc_s_fine+1)+1, ":")
                    axs[3].set_ylabel("Sat mass")
                    axs[3].set_ylim([8,12])
                    #fig.suptitle("{}: {} - {}".format(this_gal.fidx, merger.main_idx, merger.sat["idx"][0]))

                    #plt.show()
                    plt.savefig(out_base+"merger_sat_properties_{}_{}_{}.png".format(this_gal.fidx, merger.main_idx, merger.sat["idx"][0]), dpi=200)
                    #plt.close()
                    for ax in axs:
                        ax.cla()
                    axs1.cla()

       # gather data
        main_d_lambda = this_gal.finearr["lambda_r"][:-1] - this_gal.finearr["lambda_r"][1:]

        all_weights = np.zeros(len(this_gal.finearr))
        n_major_arr = np.zeros(len(this_gal.finearr))
        n_minor_arr = np.zeros(len(this_gal.finearr))

        this_gal.n_major=0
        this_gal.n_minor=0

        for merger in this_gal.mergers:
            if hasattr(merger, "i_ini"):
                for ini, end, mratio in zip(merger.i_ini, merger.i_end, merger.m_ratio):
                    #print(ini, end)
                    all_weights[end:ini+1] += mratio
                    if mratio > mratio_Major:
                        n_major_arr[end:ini+1] += 1
                        this_gal.n_major+=1
                    elif mratio > mratio_minor:
                        n_minor_arr[end:ini+1] += 1
                        this_gal.n_minor+=1

        # Merger history types.
        this_gal.n_major_arr=n_major_arr
        this_gal.n_minor_arr=n_minor_arr


        for merger in this_gal.mergers:
            if hasattr(merger, "i_ini"):
                merger.dl=[]
                #print(merger.sat["idx"][0])
                for ini, end, mratio, thisdj in zip(merger.i_ini, merger.i_end, merger.m_ratio, merger.dj):
                    merger.dl.append(np.sum(main_d_lambda[end:ini+1]))
                    merger.m_ratio_ratio = mratio/all_weights[end:ini+1]
                    dl_all.append(np.sum(main_d_lambda[end:ini+1] * (mratio/all_weights[end:ini+1])))
                    dj_all.append(thisdj)

       # for each main galaxy, compare and take larger ones, or give each merger a weight.

    if do_plot: plt.close()
    return dj_all, dl_all




def do_serialize(all_final_idxs, allresults, Allallidxs, tc,
                 nnza_cell, nnza_all,nouts,
                 prg_dir, serial_out_dir, nstep_too_short_main = 100):
    all_fid_ok=[]
    all_fidx_ok=[]

    for j, this_idx in enumerate(all_final_idxs):
        fname = prg_dir + "{}_adp.pickle".format(this_idx)
        if not os.path.isfile(fname):
            # dump_prgs broken
            continue
        print(j, "IDX=",this_idx, end="\r")
        adp = pickle.load(open(fname, "rb"))
        if min(adp[0][0]["nstep"][adp[0][0]["nstep"] > 0]) > nstep_too_short_main:
            print("skip")
            continue
        # Append age to maintree and mainresult.
        lbt = tc.zred2gyr(nnza_all.a2b(adp[0][0]["nstep"],"nstep","zred"),z_now=0)
        adp[0][0] = recf.append_fields(adp[0][0], "time", lbt)
        max_step=len(allresults)
        this_gal = serialize_results.Serial_result(adp, nnza_all)
        cnt_merger=0
        bad_main=False

        #lbt = tc.zred2gyr(nnza_all.a2b(this_gal.maintree["nstep"],"nstep","zred"), z_now=0)

        for i, this_sats in enumerate(adp):
            nout_results=[]

            for sat in this_sats:
                sat_results=[]
                for ss in sat:
                    nout=nnza_all.step2out([ss["nstep"]]) # NOT nnza_cell.
                    #print(nout)
                    if nout in nouts:
                        if i > 0:
                            if sat["nstep"] < this_gal.main_arr[-1]["nstep"]: continue
                        #print(nout, ss["idx"])
                        istep_cell = np.where(nnza_cell.nnza["nout"] == nout)[0][0]
                        allresults_now=allresults[istep_cell]
                        allresults_now_idx=Allallidxs[istep_cell]
                        i_result = np.where(allresults_now_idx == ss["idx"])[0]
                        #print("len results", len(allresults_now), "i_result", i_result)
                        if len(i_result) > 0:
                            sat_results.append(allresults_now[i_result[0]])
                            sat_results[-1].nout=int(nout)
                            sat_results[-1].nstep=nnza_cell.nnza["nstep"][istep_cell]

                #print(len(sat_results))
                nout_results.append(sat_results)
                # Merger properties
                if i == 0:
                    this_gal.main_arr = serialize_results.galresult2rec(sat_results, is_main=True)
                    #print(len(this_gal.main_arr.nstep), this_gal.main_arr.nstep.ptp())
                    if len(this_gal.main_arr.nstep) <= this_gal.main_arr.nstep.ptp():
                        #bad_main=True
                        this_gal.main_arr = serialize_results.fill_main(this_gal.main_arr, nnza_cell,tc)

                    #lbt_cell = tc.zred2gyr(nnza_cell.a2b(this_gal.main_arr["nstep"],"nstep","zred"), z_now=0)
                    this_gal.finearr = serialize_results.interpol_fine(this_gal, nnza_cell, nnza_all, tc, do_smooth=True)

                    #print("BAD", bad_main)

                elif len(sat_results) > 0 and sat_results[0].mstar > 0.0:
                    #print("merger2")

                    this_gal.add_merger(sat_results, sat)
                    cnt_merger+=1
                    #this_gal.mergers.append(serialize_results.get_merger_props(this_gal.main_arr,
                    #                        serialize_results.galresult2rec(sat_results)))
            if bad_main:
                print("Break")
                break
            #print(i)
            this_gal.data.append(nout_results)


        pickle.dump(this_gal, open(serial_out_dir+"serial_result{}.pickle".format(this_idx), "wb"))
        all_fidx_ok.append(this_idx)
        #break


    np.savetxt(serial_out_dir+"all_fidx_ok.txt", all_fidx_ok, fmt="%d")
    serial_results = [pickle.load(open(serial_out_dir+"serial_result{}.pickle".format(this_idx), "rb")) for this_idx in all_fidx_ok]

    return serial_results
