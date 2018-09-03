from utils import match as mtc
import numpy as np

def match_all_halos(gcat, hcat):
    """
    find a galaxy-halo pair.
    """
    from scipy.spatial import cKDTree

    hkdt = cKDTree(hcat.data["pos"])
    matched_halos = hmat.find_direct_halo(gcat.data, hcat.data, hkdt, n_match=100, debug=False)
    gcat.data["hosthalo"] = matched_halos["id"]

    return matched_halos


def get_kd_matches(kdtree, gal, n_match=5, rscale = 4.0, dist_upper=None):
    """
    Allow large enough rscale so that even the most sparse galaxy will find one match.

    """
    #len_tree = kdtree.length
    if dist_upper is None:
        dist_upper = gal["rvir"] * rscale

    # kdtree.query_ball_point ??
    dd, ind = kdtree.query((gal["x"],
                            gal["y"],
                            gal["z"]),
                            distance_upper_bound=dist_upper,
                            k=n_match)

    ind = ind[np.isfinite(dd)]
    dd = dd[np.isfinite(dd)]

    return dd, ind

def find_direct_halo(gdata, hdata, hkdt, n_match=50, debug=False):
    sig_pos = np.std(hdata["x"])
    sig_vel = np.std(hdata["vx"])

    matched_halo = np.zeros(len(gdata), dtype=hdata.dtype)
    matched_ids=[] # The matched_halo array is very sparse. Just keep IDs list to check for duplicates.
    miss=0
    isort = np.argsort(gdata["m"])[::-1]
    for i, thisgal in enumerate(gdata[isort]):
        r_limit = 0.8*thisgal["r"]
        dist, i_neigh = get_kd_matches(hkdt, thisgal,
                                       n_match=n_match,
                                       dist_upper=r_limit)#, dist_upper=0.5/100.)

        if debug:
            print("There are {} neighbor halos within {:.3e}".format(len(dist), r_limit))

        # Exclude already-matched haloes
        i_ok = mtc.match_list_ind(hdata["id"][i_neigh],
                                  np.setdiff1d(hdata["id"][i_neigh],
                                               matched_ids,
                                               assume_unique=True))
        dist = dist[i_ok]
        neighbor_h = hdata[i_neigh[i_ok]]

        # Mass cut.
        i_mass_ok = np.where(5*neighbor_h["mvir"] > thisgal["m"])[0]
        if debug:
            print("massive halos", np.log10(thisgal["m"]), i_mass_ok)
        neighbor_h = neighbor_h[i_mass_ok]
        dist = dist[i_mass_ok]
        if debug :
            print("Stellar mass {},  nei dist {}".format(neighbor_h["mvir"], dist))
        # 6-D dist.
        rel_vel = np.sqrt(np.einsum("...i,...i", neighbor_h["vel"] - thisgal["vel"],
                                   neighbor_h["vel"] - thisgal["vel"]))

        try:
            matched_halo[isort[i]]=neighbor_h[np.argmin(dist * rel_vel)]
            matched_ids.append(matched_halo[isort[i]]["id"])
            if debug: print("Finally.. {:.3e} \n".format(matched_halo[isort[i]]["mvir"]))
        except:
            print(i_neigh, i_ok, i_mass_ok, dist*rel_vel, "{} Failed \n".format(thisgal["id"]))
            miss+=1
            pass
    print("There are {} missing matches out of {}".format(miss, len(gdata)))

    return matched_halo
