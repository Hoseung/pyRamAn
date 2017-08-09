import utils.match as mtc
import numpy as np
from scipy.spatial import cKDTree

def find_top_host_halo(gdata, hkdt, n_match=50):
    """
    Todo
    ----
    Periodic boundary condition.
    Just copy and pad around the cube at ~ 1Mpc thickness.
    """
    matched_halo = np.zeros(len(gdata), dtype=hdata.dtype)
    for i, thisgal in enumerate(gdata[np.argsort(gdata["m"])[::-1]]):
        dist, i_neigh = get_kd_matches(hkdt, thisgal, n_match=n_match)
        touching = dist < hdata[i_neigh]["r"] + thisgal["r"]
        neighbor_h = hdata[i_neigh][touching]

        try:
            matched_halo[i]=neighbor_h[np.argmin(neighbor_h["mvir"])]
        except:
            pass

        return matched_halo

def get_kd_matches(kdtree, gal, n_match=5, rscale = 2.0, dist_upper=None):
    #len_tree = kdtree.length
    if dist_upper is None:
        dist_upper = gal["rvir"] * rscale
    dd, ind = kdtree.query((gal["x"],
                            gal["y"],
                            gal["z"]),
                            distance_upper_bound=dist_upper,
                            k=n_match)

    return dd[np.isfinite(dd)], ind[np.isfinite(dd)]


def find_direct_halo(gdata, hkdt, n_match=50):
    matched_halo = np.zeros(len(gdata), dtype=hdata.dtype)
    for i, thisgal in enumerate(gdata[np.argsort(gdata["m"])[::-1]]):
        dist, i_neigh = get_kd_matches(hkdt, thisgal, n_match=n_match)#, dist_upper=0.5/100.)

        # Exclude already-matched haloes
        id_neighbor_h_ok = np.setdiff1d(hdata["id"][i_neigh], matched_halo["id"])
        #print(id_neighbor_h_ok)
        i_ok = mtc.match_list_ind(hdata["id"][i_neigh], id_neighbor_h_ok)
        dist = dist[i_ok]
        neighbor_h = hdata[i_neigh[i_ok]]
        #print(neighbor_h["x"],neighbor_h["y"],neighbor_h["z"])

        # Mass cut.
        neighbor_h = neighbor_h[neighbor_h["m"] > thisgal["m"]]
        dist = dist[neighbor_h["m"] > thisgal["m"]]
        #print(i, len(neighbor_h), thisgal["m"], thisgal["x"], thisgal["y"], thisgal["z"])
        # 6-D dist.
        rel_vel = np.sqrt(np.einsum("...i,...i",np.column_stack((thisgal["vx"] - neighbor_h["vx"],
                                         thisgal["vy"] - neighbor_h["vy"],
                                         thisgal["vz"] - neighbor_h["vz"])),
                                        np.column_stack((thisgal["vx"] - neighbor_h["vx"],
                                         thisgal["vy"] - neighbor_h["vy"],
                                         thisgal["vz"] - neighbor_h["vz"]))))
        try:
            matched_halo[i]=neighbor_h[np.argmin(dist * rel_vel)]
        except:
            pass

    return matched_halo


def measure_density(idxs):
    """
    In comparison to Veale+17

    1) Group membership...
    2) Halo Mass
    3) Large scale density field
    4) Local galaxy density

    """
    # Halo - Galaxy Match.
    Mhal=density_halo_mass(ids, nout_fi=782, masscut = 1e10)
    #LS_d=
    Ds = density_D2N(gcat, np.array(all_fid_ok), Ns=[10,30], mass_cut=1e9)


def density_halo_mass(ids, nout_fi=782, masscut = 1e10):
    """
        IDs of halos. (Not IDxs.)

        Parameters
        ----------
        nout_fi: [=782]

        masscut:
            masscut for the halo.
            In general, no haloes smaller than the smallest galaxy are excepted.


    """
    gcat = tree.halomodule.Halo(nout=nout_fi, is_gal=True)
    hcat = tree.halomodule.Halo(nout=nout_fi, is_gal=False)


    gdata = gcat.data[mtc.match_list_ind(gcat.data["id"], ids)]
    hdata = hcat.data[np.where(hcat.data["mvir"] > masscut)[0]]
    hkdt = cKDTree(np.stack((hdata["x"], hdata["y"], hdata["z"]),axis=1))
    gkdt = cKDTree(np.stack((gdata["x"], gdata["y"], gdata["z"]),axis=1))

    halos = find_top_host_halo(gdata, hkdt, n_match=50)

    return halos["mvir"]

def density_LS_lum(gcat, rMpc=5.7):
    """
    3-D galaxy luminosity weighted density measured with a 5.7Mpc kernel.
    Well, I don't see why I should care this.
    Seems to have no merit.

    """
    raise NotImplementedError("Do I want to implement this??")
    # All galaxies.
    #gkdt = cKDTree(np.stack((gcat.data["x"], gcat.data["y"], gcat.data["z"]),axis=1))
    # Sample galaxies
    #gdata = gcat.data[mtc.match_list_ind(gcat.data["id"], IDs)]

    #DNs = np.zeros(len(gdata))
    #for i, thisgal in enumerate(gdata):
        #print(i)
    #    pass

    #return DNs


def density_D2N(gcat, IDs, Ns=[5,10], mass_cut=1e9):
    """
    Distance to the N-th nearest neighbor.

    This has been impelented elsewhere. Find it.

    Veal+17 applied Mk < −23 mag limit (MASSIVE sample).
    Cappellari+11 applied Mk < −21.5 mag limit (ATLAS sample).
    Not sure what would be the rough mass cut, but 1e4 galaxies with Mk < -23 mag are found in the MASSIVE volume (108Mpc)
    I think that's ~ 1e10.

    Now, I really have to worry about periodic boundary condition!.

    """
    n_match = max(Ns) + 1 # 0 == itself.
    # All galaxies.
    gkdt = cKDTree(np.stack((gcat.data["x"], gcat.data["y"], gcat.data["z"]),axis=1))
    # Sample galaxies
    gdata = gcat.data[mtc.match_list_ind(gcat.data["id"], IDs)]

    DNs = np.zeros((len(gdata), len(Ns)))
    for i, thisgal in enumerate(gdata):
        #print(i)
        dist, i_neigh = get_kd_matches(gkdt, thisgal, n_match=n_match, dist_upper=5./100.)
        #print(dist)
        for j, nn in enumerate(Ns):
            DNs[i,j]=dist[nn]

    return DNs
