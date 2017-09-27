import numpy as np
import pickle
import tree
from utils import hagn
from rot2.serialize_module import *
import tree.halomodule as hmo
from rot2 import cell_chunk_module as ccm

from scipy.spatial import cKDTree
from rot2.density_measure import *
from utils import cosmology
from load.info import Info
import numpy.lib.recfunctions as recf


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

def periodic_bc(data, buf=0.05):
    new=[]
    for dp1 in [1, 0, -1]:
        for dp2 in [1, 0, -1]:
            for dp3 in [1, 0, -1]:
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


def measure_P(allresults, nnza_all, nnza_cell, out_base="./RUN3/"):
    prg_dir = out_base+"test_fine_direct_prgs_gal/"
    result_dir = out_base+"lambda_results/"

    #nnza_all  = hagn.Nnza(fname = out_base + "nout_nstep_zred_aexp.txt")
    #nnza_cell = hagn.Nnza(fname = out_base + "nout_nstep_zred_aexp_63.txt")

    istep_max = 620# -500
    serial_out_dir = out_base + "result_serial/"

    #allresults = pickle.load(open(out_base+"Ids_nouts_for_P.pickle", "rb"))
    #allresults = recf.append_fields(allresults, "time", lbt)

    print("{} sample galaxies".format(len(allresults["782"])))

    info = Info(nout=787)
    tc = cosmology.Timeconvert(info, zred_now=0)
    dt = tc.zred2gyr(nnza_all.nnza["zred"][1:]) - tc.zred2gyr(nnza_all.nnza["zred"][:-1]) 

    dist_upper = 3.0 # in Mpc
    n_match=50

    nouts = nnza_all.nnza["nout"][:istep_max]

    for i, nout_now in enumerate(nouts):
        gals_now=allresults[str(nout_now)]
        info = Info(nout=nout_now)
        try:
            gdata = pickle.load(open("./GalaxyMaker/gal_pickle/gcat_{}.pickle".format(nout_now), "rb"))
            print("Good!", nout_now)
        except:
            gcat = tree.halomodule.Halo(nout=nout_now, is_gal=True)
            gdata = gcat.data

        gdata = periodic_xc(gdata)
        gkdt = cKDTree(np.stack(((gdata["x"] -0.5)*info.pboxsize,
                                 (gdata["y"] -0.5)*info.pboxsize,
                                 (gdata["z"] -0.5)*info.pboxsize), axis=1))
        P_now=[]
        #nnei =0
        for thisgal in gals_now:
            if thisgal["rgal"]==0:
                P_now.append(0)
                print("rgal==0")
                continue
            #print(thisgal["idx"])
            dist, i_neigh =  gkdt.query(gkdt.data[thisgal['id']+1],
                                         distance_upper_bound=dist_upper,
                                         k=n_match)
            i_neigh = i_neigh[np.isfinite(dist)]
            dist = dist[np.isfinite(dist)]
            #nnei += len(dist)
            #print(dist)
            # SIZE : Not real galaxy size....!! 
            neighbor_dist = dist[1:] * 1e3 # in kpc
            neighbor_mass = gdata[i_neigh[1:]]["m"]
            #print(thisgal.mstar, neighbor_mass, thisgal.rgal, neighbor_dist)
            thisgal["P_tidal"] = np.sum(neighbor_mass/thisgal["mstar"] * (thisgal["rgal"]/neighbor_dist)**3 *dt[i])
        #print(gals_now["idx"][:10])
    pickle.dump(allresults, open(out_base+"P_measured.pickle", "wb"))
    return allresults
