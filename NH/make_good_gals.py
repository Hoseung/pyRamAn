import analysis.NH_module as nhm
import tree.halomodule as hmo
import numpy as np
from analysis import NH_halo_match as hmat
import load
from load import part
from utils.sampling import Region
import utils.match as mtc
import pickle


def load_dm_of_hal(this_halo):
    reg = Region()
    reg.region_from_halo(thishalo)
    s.set_ranges(reg.ranges)

    if load_dm:
        s.add_part(ptypes=["dm id pos vel mass"])
        ind = mtc.match_list_ind(s.part.dm["id"], idlist)
        gg.dm = s.part.dm[ind]
        gg.dm["pos"] -= gg.center_code
        gg.dm["pos"] *= gg.info.boxtokpc
        gg.dm["vel"] *= gg.info.kms
        gg.dm["m"] *= gg.info.msun

def gal_hal_match(gcat, hcat, mgal_min=None, mhal_min=None):
    from scipy.spatial import cKDTree

    if mhal_min is None:
        mhal_min = mgal_min

    relevant_halos = hcat.data[hcat.data["mvir"] > mhal_min]

    kdt = cKDTree(relevant_halos["pos"])
    gals_of_interest = gcat.data[gcat.data["m"] > mgal_min]
    matched_halo = hmat.find_direct_halo(gals_of_interest, relevant_halos, kdt, n_match=100, debug=False)

    gals_of_interest["hosthalo"] = matched_halo["id"]
    gals_of_interest["host_purity"] = matched_halo["purity"]
    gals_of_interest["host_mean_m"] = matched_halo["mean_m"]

    return gals_of_interest

def cal_purity(hcat, do_part=False):
    from utils.sampling import Region
    from load import part
    import utils.match as mtc
    from scipy.spatial import cKDTree
    import pickle

    hcat.data["mean_m"] = hcat.data["m"] / hcat.data["np"]
    min_mean_m = hcat.data["mean_m"].min()

    if do_part:
        reg = Region()

        # Part alone should work!!
        ptcl = part.Part(info=hcat.info, ptypes=["dm id pos mass"])
        h_center = hcat.data[np.argmax(hcat.data["np"])]
        reg.region_from_halo(h_center)
        reg.radius = 0.1
        ptcl.set_ranges(reg.ranges)
        ptcl.load()
        pkdt = cKDTree(ptcl.dm["pos"])

        M_hires = 1.3e6 / ptcl.info.msun

        for i, (this_halo, this_idlist) in enumerate(zip(hcat.data, hcat.idlists)):

            i_near = pkdt.query_ball_point((this_halo["x"], this_halo["y"], this_halo["z"]),
                                            this_halo["r"] * 1.5)
            dd = ptcl.dm[i_near]
            dm = dd[mtc.match_list_ind(dd["id"], this_idlist, allow_swap=False)]
            # print("contam: {:.2f} %".format(100 * np.sum(dm["m"] > M_hires)/ len(dm)))
            this_halo["purity"] = 1 - np.sum(dm["m"] > M_hires) / len(dm)
            this_halo["np"] = len(dm)
    else:
        hcat.data["purity"] = 2 - hcat.data["mean_m"]/min_mean_m

    output = np.vstack((hcat.data["id"], hcat.data["purity"],
                        hcat.data["mean_m"], hcat.data["np"])).T
    pickle.dump(output, open("purity_{}.pickle".format(hcat.nout), "wb"))
    np.savetxt("purity_{}.txt".format(hcat.nout), output, fmt="  %5d   %.5f   %.5e    %7d")

    #return purity


def make_match_gal_hal_zoom(gcat, hcat, purity_crit=0.99):
    """
    Match galaxy and "pure" halos.

    parameters
    ----------
    purity_crit : {0.99}
        Contamination of a halo must be less than 1-purity_crit.
    """
    # temporarily modify gcat fields
    assert gcat.nout == hcat.nout, "galxy catalog and halo catalog are from different nout"
    nout = hcat.nout
    lnames = list(gcat.data.dtype.names)
    lnames[lnames.index("sig")] = "host_purity"
    lnames[lnames.index("sigbulge")] = "host_mean_m"
    gcat.data.dtype.names = tuple(lnames)

    # temporarily modify hcat fields
    lnames = list(hcat.data.dtype.names)
    lnames[lnames.index("p_rho")] = "purity"
    lnames[lnames.index("p_c")] = "mean_m"
    hcat.data.dtype.names = tuple(lnames)

    cal_purity(hcat)
    gals_of_interest = gal_hal_match(gcat, hcat, mgal_min=3e8)
    good_gals = gals_of_interest[gals_of_interest["host_purity"] > purity_crit]

    pickle.dump(good_gals, open("good_gals_{}.pickle".format(nout), "wb"))
    # release memory... how?
    gcat = None


if __name__=='__main__':
    out_base='./'
    wdir = "./"

    nouts = np.arange(594, 593, -1)
    for nout in nouts:
        print(" \n NOUT = {}\n".format(nout))
        # NH galaxy - Double
        # NH halo   - Single
        gcat = hmo.Halo(nout=nout, is_gal=True, double=True)
        hcat = hmo.Halo(nout=nout, is_gal=False, double=False, return_id=True)
        make_match_gal_hal_zoom(gcat, hcat, purity_crit=0.995)
