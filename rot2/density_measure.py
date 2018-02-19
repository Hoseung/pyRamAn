import utils.match as mtc
import numpy as np
from scipy.spatial import cKDTree

def periodic_bc(data, buf=0.05):
    """
    should run for 0,0,0 first so that index and id match remains.
    """
    new=[]
    for dp1 in [0, 1, -1]:
        for dp2 in [0,1, -1]:
            for dp3 in [0,1, -1]:
                temp = data.copy()
                if dp1 != 0:
                    temp["x"]+=dp1
                if dp2 != 0:
                    temp["y"]+=dp2
                if dp3 != 0:
                    temp["z"]+=dp3
                new.append(temp)

    new = np.concatenate(new)
    return new[np.where( (new["x"] > -buf)*(new["x"] < 1+buf)*
                         (new["y"] > -buf)*(new["y"] < 1+buf)*
                         (new["z"] > -buf)*(new["z"] < 1+buf))[0]]


def get_kd_matches(kdtree, gal, n_match=5, rscale = 4.0, dist_upper=None):
    """
    Allow large enough rscale so that even the most sparse galaxy will find one match.

    """
    #len_tree = kdtree.length
    if dist_upper is None:
        dist_upper = gal["rvir"] * rscale
    dd, ind = kdtree.query((gal["x"],
                            gal["y"],
                            gal["z"]),
                            distance_upper_bound=dist_upper,
                            k=n_match)

    ind = ind[np.isfinite(dd)]
    dd = dd[np.isfinite(dd)]

    return dd,ind


def find_direct_halo(gdata, hdata, hkdt, n_match=50):
    matched_halo = np.zeros(len(gdata), dtype=hdata.dtype)
    miss=0
    isort = np.argsort(gdata["m"])[::-1]
    for i, thisgal in enumerate(gdata[isort]):
        dist, i_neigh = get_kd_matches(hkdt, thisgal, n_match=n_match)#, dist_upper=0.5/100.)

        # Exclude already-matched haloes
        id_neighbor_h_ok = np.setdiff1d(hdata["id"][i_neigh], matched_halo["id"])

        i_ok = mtc.match_list_ind(hdata["id"][i_neigh], id_neighbor_h_ok)
        dist = dist[i_ok]
        neighbor_h = hdata[i_neigh[i_ok]]

        # Mass cut.
        i_mass_ok =neighbor_h["m"] > thisgal["m"]
        neighbor_h = neighbor_h[i_mass_ok]
        dist = dist[i_mass_ok]
        # 6-D dist.
        rel_vel = np.sqrt(np.einsum("...i,...i",np.column_stack((thisgal["vx"] - neighbor_h["vx"],
                                         thisgal["vy"] - neighbor_h["vy"],
                                         thisgal["vz"] - neighbor_h["vz"])),
                                        np.column_stack((thisgal["vx"] - neighbor_h["vx"],
                                         thisgal["vy"] - neighbor_h["vy"],
                                         thisgal["vz"] - neighbor_h["vz"]))))
        try:
            matched_halo[isort[i]]=neighbor_h[np.argmin(dist * rel_vel)]
        except:
            miss+=1
            pass
    print("There are {} missing matches".format(miss))

    return matched_halo

def find_top_host_halo(gdata, hdata, hkdt, n_match=50, rscale=1.0):
    """
    Todo
    ----
    Periodic boundary condition.
    Just copy and pad around the cube at ~ 1Mpc thickness.
    """
    matched_halo = np.zeros(len(gdata), dtype=hdata.dtype)
    isort = np.argsort(gdata["m"])[::-1]
    for i, thisgal in enumerate(gdata[isort]):
        dist, i_neigh = get_kd_matches(hkdt, thisgal, n_match=n_match, dist_upper=0.2)
        touching = dist < (hdata[i_neigh]["r"] + thisgal["r"]) * rscale
        neighbor_h = hdata[i_neigh[touching]]
        #print(dist[touching])

        try:
        #if True:
            matched_halo[isort[i]]=neighbor_h[np.argmax(neighbor_h["mvir"])]
        except:
            print("No Most massive halo...???")
            pass

    return matched_halo


def match_halo_gal(ids, gcdata, hcdata, masscut=None):
    """
        IDs of halos. (Not IDxs.)

        Parameters
        ----------

        masscut:
            masscut for the halo.
            In general, no haloes smaller than the smallest galaxy are excepted.

    """
    gdata = gcdata[mtc.match_list_ind(gcdata["id"], ids)]
    mcut = np.median(gdata["m"]) * 0.2 # Host halo must be larger than 20% of Mstar.
    hdata = hcdata[np.where(hcdata["mvir"] > mcut)[0]]
    hdata = periodic_bc(hdata)
    hkdt = cKDTree(np.stack((hdata["x"], hdata["y"], hdata["z"]),axis=1))
    gkdt = cKDTree(np.stack((gdata["x"], gdata["y"], gdata["z"]),axis=1))

    direct_halos = find_direct_halo(gdata, hdata, hkdt, n_match=500)
    halos = find_top_host_halo(gdata, hdata, hkdt, n_match=2000)
    halos2 = find_top_host_halo(gdata, hdata, hkdt, n_match=3000, rscale=2.0)

    return direct_halos, halos, halos2


def measure_density(idxs):
    """
    In comparison to Veale+17

    1) Group membership...
    2) Halo Mass
    3) Large scale density field
    4) Local galaxy density

    """
    # Halo - Galaxy Match.
    Mhal=density_halo_mass(ids, masscut = 1e10)
    #LS_d=
    Ds = density_D2N(gcat, np.array(all_fid_ok), Ns=[10,30])


def density_halo_mass(ids, gcat, hcat, masscut = 1e10):
    """
        IDs of halos. (Not IDxs.)

        Parameters
        ----------

        masscut:
            masscut for the halo.
            In general, no haloes smaller than the smallest galaxy are excepted.


    """
    gdata = gcat.data[mtc.match_list_ind(gcat.data["id"], ids)]
    hdata = hcat.data[np.where(hcat.data["mvir"] > masscut)[0]]
    hkdt = cKDTree(np.stack((hdata["x"], hdata["y"], hdata["z"]),axis=1))

    # !!!!!! find_direct_halo sorts galaxy order. check for consistency.
    direct_halos = find_direct_halo(gdata, hdata, hkdt, n_match=1000)
    halos = find_top_host_halo(gdata, hdata, hkdt, n_match=1000)

    return direct_halos["mvir"], halos["mvir"]


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

def density_D2N(gcdata, info, serial_results, Ns=[10,50], dist_upper=25.):
    """
    Distance to the N-th nearest neighbor.

    This has been impelented elsewhere. Find it.

    Veal+17 applied Mk < −23 mag limit (MASSIVE sample).
    Cappellari+11 applied Mk < −21.5 mag limit (ATLAS sample).
    Not sure what would be the rough mass cut, but 1e4 galaxies with Mk < -23 mag are found in the MASSIVE volume (108Mpc)
    I think that's ~ 1e10.

    Now, I really have to worry about periodic boundary condition!.

    """
    nout = info.nout
    pb = info.pboxsize

    n_match = 2 * max(Ns) + 1 # 0 == itself.

    # All galaxies.
    gdata = gcdata[gcdata["m"] >= 1e10]
    gdata = periodic_bc(gdata)
    gkdt = cKDTree(np.stack(((gdata["x"]-0.5)*pb,
                             (gdata["y"]-0.5)*pb,
                             (gdata["z"]-0.5)*pb),axis=1))
    # Sample galaxies
    for tg in serial_results:
        vals = tg.finedata
        inow=np.where(vals["nout"]==nout)[0]
        if len(inow) > 0:
            val_now = vals[inow][0]
            dist, i_neigh =  gkdt.query(val_now["pos"],
                                 distance_upper_bound=dist_upper,
                                 k=n_match)
#            print("gal", val_now["pos"])
            i_neigh = i_neigh[np.isfinite(dist)]
            dist = dist[np.isfinite(dist)]
            for j, nn in enumerate(Ns):
                tg.env["d5"][inow]=dist[5]
                tg.env["d10"][inow]=dist[10]
                tg.env["d50"][inow]=dist[50]


def find_top_host_halo_old(gdata, hdata, hkdt, n_match=50, rscale=1.0):
    """
    Todo
    ----
    Periodic boundary condition.
    Just copy and pad around the cube at ~ 1Mpc thickness.
    """
    matched_halo = np.zeros(len(gdata), dtype=hdata.dtype)
    for i, thisgal in enumerate(gdata[np.argsort(gdata["m"])[::-1]]):
        dist, i_neigh = get_kd_matches(hkdt, thisgal, n_match=n_match, dist_upper=0.1)
        touching = dist < (hdata[i_neigh]["r"] + thisgal["r"]) * rscale
        neighbor_h = hdata[i_neigh][touching]

        try:
        #if True:
            matched_halo[i]=neighbor_h[np.argmax(neighbor_h["mvir"])]
        except:
            print("No Most massive halo...???")
            pass

    return matched_halo

def density_D2N_old(gcdata, IDs, Ns=[5,10], dist_upper=15./100.):
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
    gdata_target = gcdata[mtc.match_list_ind(gcdata["id"], IDs)]

    # All galaxies.
    #print("len_before", len(gcat.data))
    gdata = gcdata[gcdata["m"] >= min(gdata_target["m"])]
    gdata = periodic_bc(gdata)
    #print("len_after", len(gdata))
    gkdt = cKDTree(np.stack((gdata["x"], gdata["y"], gdata["z"]),axis=1))
    # Sample galaxies

    DNs = np.zeros((len(gdata_target), len(Ns)))
    for i, thisgal in enumerate(gdata_target):
        dist, i_neigh = get_kd_matches(gkdt, thisgal, n_match=n_match, dist_upper=dist_upper)
        for j, nn in enumerate(Ns):
            DNs[i,j]=dist[nn]

    return DNs


def measure_P(gdata, info, serial_results, dt,
                                dist_upper = 3.0,
                                n_match=1000, short=False):
    """

        Dist upper in Mpc.
    """
    nout = info.nout

    #gdata = gcat.data[mtc.match_list_ind(gcat.data["id"], ids)]
    gdata = periodic_bc(gdata)

    gkdt = cKDTree(np.stack(((gdata["x"] -0.5)*info.pboxsize,
                             (gdata["y"] -0.5)*info.pboxsize,
                             (gdata["z"] -0.5)*info.pboxsize), axis=1))

    for thisgal in serial_results:
        if short:
            vals = thisgal.main_data
        else:
            vals = thisgal.finedata
        inow=np.where(vals["nout"]==nout)[0]
        if len(inow) > 0:
            if vals["rgal"][inow]==0:
                print("rgal==0")
                continue
            else:
                vals_now = vals[inow]
                dist, i_neigh =  gkdt.query(vals_now['pos'],
                                     distance_upper_bound=dist_upper,
                                     k=n_match)
                i_neigh = i_neigh[np.isfinite(dist)]
                dist = dist[np.isfinite(dist)]
                neighbor_dist = dist[1:] * 1e3 # in kpc
                neighbor_mass = gdata[i_neigh[1:]]["m"]
                #thisgal["P_tidal_h_nout"] = nout_now
                if short:
                    thisgal.env_short["P_tidal_h"][inow] = np.sum(neighbor_mass/vals_now["mstar"] \
                                                   * (vals_now["rgal"]/neighbor_dist)**3) *dt
                else:
                    thisgal.env["P_tidal"][inow] = np.sum(neighbor_mass/vals_now["mstar"] \
                                                   * (vals_now["rgal"]/neighbor_dist)**3) *dt
        else:
            print("No matching galaxy at {}".format(nout))
