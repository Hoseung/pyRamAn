import numpy as np
import pickle
from utils import hagn
from rot2.analysis import *
import tree.halomodule as hmo
from rot2 import cell_chunk_module as ccm
from scipy.spatial import cKDTree
from rot2.density_measure import *
from utils import cosmology
from load.info import Info


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

if __name__ == "__main__":
    out_base="./RUN1/"

    prg_dir = out_base+"test_direct_prgs_gal/"
    result_dir = out_base+"lambda_results/"

    # Load prg tree data
    all_ids= np.genfromtxt(prg_dir + "final_idxs_allmassive_gal.txt", dtype=int)
    all_final_idxs = all_ids[:,0]
    all_final_ids = all_ids[:,1]

    nnza_all = hagn.Nnza(fname=out_base+"nout_nstep_zred_aexp.txt")

    istep_max = 700
    nouts = nnza_all.nnza["nout"][:istep_max]

    info = Info(nout=787)
    tc = cosmology.Timeconvert(info, zred_now=0)

    dt = tc.zred2gyr(nnza_all.nnza["zred"][1:]) - tc.zred2gyr(nnza_all.nnza["zred"][:-1])

    Ps = []

    dist_upper = 3.0 # in Mpc
    n_match=50
    xs=[]
    ys=[]
    zs=[]

    ngal_all = 0
    for this_gal in serial_results:
        ngal_all+=len(this_gal.finearr)

    all_gal_prps = np.zeros(ngal_all, dtype=[("nout", "<i4"),
                                             ("id", "<i8"),
                                             ("mass", "<f8"),
                                             ("rgal", "<f8"),
                                             ("pos", "<f8", (3)),
                                             ("P", "<f8")])
    ind=0
    for this_gal in serial_results:
        ll = len(this_gal.finearr)
        all_gal_prps[ind:ind+ll]["nout"] = this_gal.finearr["nout"]
        all_gal_prps[ind:ind+ll]["id"] = this_gal.finearr["id"]
        all_gal_prps[ind:ind+ll]["mass"] = this_gal.finearr["mass"]
        all_gal_prps[ind:ind+ll]["rgal"] = this_gal.finearr["rgal"]
        all_gal_prps[ind:ind+ll]["pos"] = this_gal.finearr["pos"]
        ind +=ll

    iall=0
    for i, nout_now in enumerate(np.unique(all_gal_prps["nout"])):
        i_now = np.where(all_gal_prps["nout"] == nout_now)[0]

        gals_now = all_gal_prps[i_now]
        #nout_now = nouts[i]
        nout_now = int(nout_now)
        info = Info(nout=nout_now)
        gcat = tree.halomodule.Halo(nout=nout_now, is_gal=True)
        gdata = gcat.data
        #gdata = pickle.load(open("./GalaxyMaker/gal/gcat_{}.pickle".format(nout_now), "rb"))
        gkdt = cKDTree(np.stack(((gdata["x"] -0.5)*info.pboxsize,
                                 (gdata["y"] -0.5)*info.pboxsize,
                                 (gdata["z"] -0.5)*info.pboxsize), axis=1))

        P_now=[]
        for thisgal in gals_now:
            if thisgal["rgal"] == 0:
                P_now.append(0)
                continue

            dist, i_neigh =  gkdt.query( thisgal["pos"],
                                         distance_upper_bound=dist_upper,
                                         k=n_match)
            i_neigh = i_neigh[np.isfinite(dist)]
            dist = dist[np.isfinite(dist)]

            # SIZE : Not real galaxy size....!!
            neighbor_dist = dist[1:] * 1e3 # in kpc
            neighbor_mass = gdata[i_neigh[1:]]["m"]
            #print(thisgal.mstar, neighbor_mass, thisgal.rgal, neighbor_dist)

            #Ps.append()[iall] = (nout_now,
            #         thisgal_id,
            #         np.sum(thisgal.mstar/neighbor_mass * (thisgal.rgal/neighbor_dist)**3 *dt[i]))
            P_now.append(np.sum(thisgal["mass"]/neighbor_mass * (thisgal["rgal"]/neighbor_dist)**3 *dt[i]))
            #iall+=1
        #Ps.append(P_now)
        all_gal_prps[i_now] = P_now

    pickle.dump(all_gal_prps, open("P_measured.pickle", "wb"))
