import numpy as np
import utils.match as mtc
from scipy.signal import savgol_filter
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d, Akima1DInterpolator
import matplotlib.pyplot as plt

def add_main(md, data, ind):
    md[ind]["mstar"] = data.mstar
    md[ind]["pos"] = (data.xc, data.yc, data.zc)
    md[ind]["vel"] =(data.vxc, data.vyc, data.vzc)
    if hasattr(data, "mge_result_list"):
        md[ind]["eps"] = data.mge_result_list[0]["eps"]
    try:
        md[ind]["lambda_r"] = data.lambda_result_list[0][4]
    except:
        md[ind]["lambda_r"] = data.lambda_r

    if data.lvec is not None:
        md[ind]["reff"] = data.reff
        md[ind]["nvec"] = data.nvec
        # If there is Lvec, there is rgal.
        md[ind]["lvec"] = data.lvec
        md[ind]["rgal"] = data.rgal

    if hasattr(data,"gas_results"):
        md[ind]["mgas"] = data.gas_results["mgas_tot"]
        md[ind]["mgas_cold"] = data.gas_results["mgas_cold"]
        md[ind]["lgas"]= data.gas_results["Ln_gas"]

    try:
        md[ind]["vsig"]=data.vsig_results["V_sig"]
    except:
        md[ind]["vsig"]=0

    if hasattr(data, "sfr_results"):
        md[ind]["sfr01"]=data.sfr_results["sfrs"][0]
        md[ind]["sfr05"]=data.sfr_results["sfrs"][1]
        md[ind]["sfr1"]=data.sfr_results["sfrs"][2]
        md[ind]["sfr_area"]=data.sfr_results["area"]


# Interpolation tools
def interp(x1, y1, x2, too_many_nan_frac=0.3):
    """
    often angle measurements have nan.
    Give zero weight to nan values.

    Giving zero weight to nans works with nans in the array.
    But cases with leading or trailing nans are not solved.
    Just opt out nans.
    """
    w = np.isnan(y1).astype(bool)# * (y1 == 0.0)
    if np.sum(w) > 0.3* len(x1):
        raise ValueError("Too many nans in the array during interpolation")
    f = InterpolatedUnivariateSpline(x1[~w],y1[~w])#, w=~w[::-1])
    return f(x2)


def interp_akima(x1, y1, x2, too_many_nan_frac=0.3):
    w = np.isnan(y1).astype(bool)# * (y1 == 0.0)
    if np.sum(w) > 0.3* len(x1):
        raise ValueError("Too many nans in the array during interpolation")
    f = Akima1DInterpolator(x1[~w], y1[~w])
    return f(x2)


def interp_inv(x1, y1, x2, too_many_nan_frac=0.3):
    """
    often angle measurements have nan.
    Give zero weight to nan values.

    Giving zero weight to nans works with nans in the array.
    But cases with leading or trailing nans are not solved.
    Just opt out nans.
    """
    w = np.isnan(y1).astype(bool)# * (y1 == 0.0)
    if np.sum(w) > 0.3* len(x1):
        raise ValueError("Too many nans in the array during interpolation")
    #print(y1[~w][::-1])
    f = InterpolatedUnivariateSpline(x1[~w][::-1],y1[~w][::-1])#, w=~w[::-1])
    return f(x2)


def interp_np(x1, y1, x2, too_many_nan_frac=0.4):
    """
    called when the length is <= 5.
    So, 2 out of 5 is 40% of nan.

    """
    w = np.isnan(y1).astype(bool)
    if np.sum(w) > too_many_nan_frac * len(x1):
        raise ValueError("Too many nans in the array during interpolation")

    return np.interp(x2, x1[~w], y1[~w])

def interp_np_inv(x1, y1, x2, too_many_nan_frac=0.4):
    """
    called when the length is <= 5.
    So, 2 out of 5 is 40% of nan.

    """
    w = np.isnan(y1).astype(bool)
    if np.sum(w) > too_many_nan_frac * len(x1):
        raise ValueError("Too many nans in the array during interpolation")

    return np.interp(x2[::-1], x1[~w][::-1], y1[~w][::-1])[::-1]



def norm_3d(vec):
    return vec/np.sqrt(np.einsum("...i,...i", vec, vec))[:,None]

def interp_vec_arr_to_arrays(nout, vec, new_nout, normalize=True):
    """
    normal vectors must have the magnitude 1.
    but interpolation can result in values slightly larger than 1.
    So, normalized it again.
    """
    #print("before interpol", vec[:,0])
    if len(vec) > 5:
        vecx = interp(nout, vec[:,0], new_nout)
        vecy = interp(nout, vec[:,1], new_nout)
        vecz = interp(nout, vec[:,2], new_nout)
    else:
        vecx = interp_np(nout, vec[:,0], new_nout)
        vecy = interp_np(nout, vec[:,1], new_nout)
        vecz = interp_np(nout, vec[:,2], new_nout)
        #vecx = np.interp(new_nout[::-1], nout[::-1], vec[::-1,0])[::-1]
        #vecy = np.interp(new_nout[::-1], nout[::-1], vec[::-1,1])[::-1]
        #vecz = np.interp(new_nout[::-1], nout[::-1], vec[::-1,2])[::-1]
    #print("After interpol", vecx)
    if normalize:
        mag = 1./(np.sqrt(vecx**2+vecy**2+vecz**2))
        vecx *= mag
        vecy *= mag
        vecz *= mag

    return vecx,vecy,vecz


def smooth2(x, window_length, polyorder):
    savgol_filter(x, window_length, polyorder, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)

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

def get_sat_data_fine(sat_data ):
    """
    data, nout, new_nout,,, how do I know new_nout...??
    """

    pass

def get_main_data_fine(tg,
                       win_size=15,
                       pol_deg=2):

    md = tg.main_data
    # require mstar > 1e6 to be a detection.
    #print(md["mstar"][::-1])
    i_non_zero = len(md) - np.argmax(md["mstar"][::-1] > 1e6) -1
    md_lbt_max = md["lbt"][i_non_zero]
    mt = tg.finedata[tg.finedata["lbt"] <= md_lbt_max]

    newarr = np.zeros(len(mt), dtype=md.dtype)
    newarr["nout"] = mt["nout"]
    newarr["lbt"]  = mt["lbt"]
    newarr["mstar"] = smooth(interp_np(md["lbt"], md["mstar"], mt["lbt"]))
    newarr["mgas"] = smooth(interp_np(md["lbt"], md["mgas"], mt["lbt"]))
    newarr["mgas_cold"] = smooth(interp_np(md["lbt"], md["mgas_cold"], mt["lbt"]))
    lx,ly,lz = interp_vec_arr_to_arrays(md["lbt"], md["lgas"], mt["lbt"])
    newarr["lgas"][:,0] = savgol_filter(lx, win_size, pol_deg)
    newarr["lgas"][:,1] = savgol_filter(ly, win_size, pol_deg)
    newarr["lgas"][:,2] = savgol_filter(lz, win_size, pol_deg)

    newarr["rgal"] = smooth(interp_np(md["lbt"], md["rgal"], mt["lbt"]))
    newarr["reff"] = smooth(interp_np(md["lbt"], md["reff"], mt["lbt"]))

    newarr["sfr01"] = savgol_filter(interp(md["lbt"], md["sfr01"], mt["lbt"]),
                                    win_size, pol_deg)
    newarr["sfr05"] = savgol_filter(interp(md["lbt"], md["sfr05"], mt["lbt"]),
                                    win_size, pol_deg)
    newarr["sfr1"]  = savgol_filter(interp(md["lbt"],  md["sfr1"], mt["lbt"]),
                                    win_size, pol_deg)
    newarr["sfr_area"] = savgol_filter(interp(md["lbt"],  md["sfr_area"], mt["lbt"]),
                                    win_size, pol_deg)
    print("interp_akima")
    newarr["eps"] = interp_akima(md["lbt"],  md["eps"], mt["lbt"])
                    #savgol_filter(interp(md["lbt"],  md["eps"], mt["lbt"]),
                    #                win_size, pol_deg)
    newarr["lambda_r"] = interp_akima(md["lbt"],  md["lambda_r"], mt["lbt"])
        #savgol_filter(interp(md["lbt"],  md["lambda_r"], mt["lbt"]),
                         #           win_size, pol_deg)
    newarr["vsig"]  = savgol_filter(interp(md["lbt"],  md["vsig"], mt["lbt"]),
                                    win_size, pol_deg)

    lx,ly,lz = interp_vec_arr_to_arrays(md["lbt"], md["lvec"], mt["lbt"])
    newarr["lvec"][:,0] = savgol_filter(lx, win_size, pol_deg)
    newarr["lvec"][:,1] = savgol_filter(ly, win_size, pol_deg)
    newarr["lvec"][:,2] = savgol_filter(lz, win_size, pol_deg)
    newarr["nvec"] = norm_3d(newarr["lvec"])

    return newarr


def main_data_to_tree(tg,
                      win_size=15,
                      pol_deg=2,
                      akima=False):

    md = tg.main_data
    mt = tg.finedata
    # require mstar > 1e6 to be a detection.
    md_lbt_max = md["lbt"][len(md) - np.argmax(md["mstar"][::-1] > 1e6) -1]
    ind = np.where(mt["lbt"] <= md_lbt_max)[0]
    lbt_new=mt["lbt"][ind]

    mt["mstar"][ind] = smooth(interp_np(md["lbt"], md["mstar"], lbt_new))
    mt["mgas"][ind] = smooth(interp_np(md["lbt"], md["mgas"], lbt_new))
    mt["mgas_cold"][ind] = smooth(interp_np(md["lbt"], md["mgas_cold"], lbt_new))
    lx,ly,lz = interp_vec_arr_to_arrays(md["lbt"], md["lgas"], lbt_new)
    mt["lgas"][ind,0] = savgol_filter(lx, win_size, pol_deg)
    mt["lgas"][ind,1] = savgol_filter(ly, win_size, pol_deg)
    mt["lgas"][ind,2] = savgol_filter(lz, win_size, pol_deg)

    mt["rgal"][ind] = smooth(interp_np(md["lbt"], md["rgal"], lbt_new))
    mt["reff"][ind] = smooth(interp_np(md["lbt"], md["reff"], lbt_new))

    mt["sfr01"][ind] = savgol_filter(interp(md["lbt"], md["sfr01"], lbt_new),
                                    win_size, pol_deg)
    mt["sfr05"][ind] = savgol_filter(interp(md["lbt"], md["sfr05"], lbt_new),
                                    win_size, pol_deg)
    mt["sfr1"][ind]  = savgol_filter(interp(md["lbt"],  md["sfr1"], lbt_new),
                                    win_size, pol_deg)
    mt["sfr_area"][ind] = savgol_filter(interp(md["lbt"],  md["sfr_area"], lbt_new),
                                    win_size, pol_deg)
    if akima:
        mt["eps"][ind]  = interp_akima(md["lbt"],  md["eps"], lbt_new)
        mt["lambda_r"][ind] = interp_akima(md["lbt"],  md["lambda_r"], lbt_new)
    else:
        mt["eps"][ind]  = savgol_filter(interp(md["lbt"],  md["eps"], lbt_new),
                                        win_size, pol_deg)
        mt["lambda_r"][ind] = savgol_filter(interp(md["lbt"],  md["lambda_r"], lbt_new),
                                        win_size, pol_deg)
    mt["vsig"][ind] = savgol_filter(interp(md["lbt"],  md["vsig"], lbt_new),
                                    win_size, pol_deg)

    lx,ly,lz = interp_vec_arr_to_arrays(md["lbt"], md["lvec"], lbt_new)
    mt["lvec"][ind,0] = savgol_filter(lx, win_size, pol_deg)
    mt["lvec"][ind,1] = savgol_filter(ly, win_size, pol_deg)
    mt["lvec"][ind,2] = savgol_filter(lz, win_size, pol_deg)
    mt["nvec"][ind] = norm_3d(mt["lvec"][ind])


def check_main_fine(data, fidx, suffix=""):
    fig, axs = plt.subplots(3,4)
    fig.set_size_inches(16,12)

    axs[0,0].plot(data["nout"], np.log10(data["mstar"]))
    axs[0,0].set_title("log mass")
    axs[0,1].plot(data["nout"], np.log10(data["mgas"] + 1e4))
    axs[0,1].set_title("mgas")
    axs[0,2].plot(data["nout"], np.log10(data["mgas_cold"] + 1e4))
    axs[0,2].set_title("mgas_cold")
    axs[0,3].plot(data["nout"], data["lgas"][:,0], label="lx")
    axs[0,3].plot(data["nout"], data["lgas"][:,1], label="ly")
    axs[0,3].plot(data["nout"], data["lgas"][:,2], label="lz")
    axs[0,3].set_title("gas Lvec")

    axs[1,0].plot(data["nout"], data["rgal"])
    axs[1,0].set_title("Rgal")
    axs[1,1].plot(data["nout"], data["lambda_r"])
    axs[1,1].set_title("lambda")
    axs[1,2].plot(data["nout"], data["eps"])
    axs[1,2].set_title("eps")
    axs[1,3].plot(data["nout"], data["vsig"])
    axs[1,3].set_title("V/sig")

    axs[2,0].plot(data["nout"], data["sfr01"], label="sfr01")
    axs[2,0].plot(data["nout"], data["sfr05"], label="sfr05")
    axs[2,0].plot(data["nout"], data["sfr1"], label="sfr1")
    axs[2,0].set_title("sfr")
    axs[2,1].plot(data["nout"], data["lvec"][:,0], label="lx")
    axs[2,1].plot(data["nout"], data["lvec"][:,1], label="ly")
    axs[2,1].plot(data["nout"], data["lvec"][:,2], label="lz")
    axs[2,1].set_title("Stellar Lvec")
    axs[2,2].plot(data["nout"], data["P_tidal_h"])
    axs[2,2].set_title("P_tidal_h")
    axs[2,3].plot(data["nout"], data["P_tidal_h"])
    axs[2,3].set_title("P_tidal_h")

    plt.savefig("{}_{}.png".format(fidx, suffix))
    plt.close()


# check al variables
def check_merger_prps(merger, hidx, midx, suffix=""):
    fig, axs = plt.subplots(3,4)
    fig.set_size_inches(16,12)
    axs[0,0].plot(merger["nout"], merger["dist"])
    axs[0,0].set_title("dist")
    axs[0,1].plot(merger["nout"], np.log10(merger["mstar"]))
    axs[0,1].set_title("log mass")
    axs[0,2].plot(merger["nout"], merger["m_frac"])
    axs[0,2].set_title("m_frac")

    axs[1,0].plot(merger["nout"], merger["rgal_sum"])
    axs[1,0].set_title("Rgal sum")
    axs[1,1].plot(merger["nout"], np.sqrt(np.einsum("...i,...i", merger["jorbit"], merger["jorbit"])))
    axs[1,1].set_title("jorbit mag")
    axs[1,2].plot(merger["nout"], merger["orbitang"])
    axs[1,2].set_title("orbitang")
    axs[1,3].plot(merger["nout"], merger["spinang"])
    axs[1,3].set_title("spinang")

    axs[2,0].plot(merger["nout"], merger["spinmag"])
    axs[2,0].set_title("spin mag")
    axs[2,1].plot(merger["nout"], np.sqrt(np.einsum("...i,...i", merger["rel_vel"], merger["rel_vel"])))
    axs[2,1].set_title("rel_vel mag")
    axs[2,2].plot(merger["nout"], merger["reff_s"])
    axs[2,2].set_title("reff")

    plt.savefig("{}_{}{}.png".format(hidx, midx, suffix))
    plt.close()
