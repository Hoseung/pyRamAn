import analysis.NH_module as nhm
import tree.halomodule as hmo 
import numpy as np
from analysis import NH_halo_match as hmat
import load
from load import part
from utils.sampling import Region

import importlib
import utils.match as mtc


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

def cal_purity(hcat):
    from utils.sampling import Region
    from load import part
    import utils.match as mtc

    reg = Region()

    # Part alone should work!!
    ptcl = part.Part(info=hcat.info, ptypes=["dm id pos mass"])
    M_hires = 1.3e6 / ptcl.info.msun

    hcat.data["mean_m"] = hcat.data["m"] / hcat.data["np"]
    min_mean_m = hcat.data["mean_m"].min()
    i_ok = np.where(hcat.data["mean_m"] < 2*min_mean_m)[0]

    for i, (this_halo, this_idlist) in enumerate(zip(hcat.data, hcat.idlists)):
        dist = np.sqrt(np.sum(np.square(this_halo["pos"] - [0.49683937, 0.50610047, 0.51450588])))
        if dist > 0.2:
            this_halo["purity"] = 0
            continue
        print("Dist from center", dist)
        
        reg.region_from_halo(this_halo)
        reg.radius *= 1.2
        ptcl.set_ranges(reg.ranges)

        ptcl.load(return_whole=False)
        
        dm = ptcl.dm[mtc.match_list_ind(ptcl.dm["id"], this_idlist, allow_swap=False)]
        print("contam:", np.sum(dm["m"] > M_hires), len(dm))
        this_halo["purity"] = 1 - np.sum(dm["m"] > M_hires) / len(dm)
        this_halo["np"] = len(dm)
        
    np.savetxt("purity_{}.txt".format(nout),
         np.vstack((hcat.data["id"], hcat.data["purity"], hcat.data["mean_m"], hcat.data["np"]),
                    fmt="  %5s   %.5f   %.5f    %7d"))
    #pickle.dump(purity, open("purity_{}.pickle".format(nout), "wb"))

    #return purity


if __name__=='__main__':
    out_base='./'
    wdir = "./"

    nouts = np.arange(629, 29, -1)

    for nout in nouts:    
        gcat = hmo.Halo(nout=nout, is_gal=True, double=True, pure=False)
        hcat = hmo.Halo(nout=nout, is_gal=False, double=False, pure=False, return_id=True)
        
        # temporarily modify gcat fields 
        lnames = list(gcat.data.dtypes.names)
        lnames[lnames.index("sig")] = "host_purity"
        lnames[lnames.index("sigbulge")] = "host_mean_m"
        gcat.data.dtypes = tuple(lnames)
        
        # temporarily modify hcat fields
        lnames = list(hcat.data.dtypes.names)
        lnames[lnames.index("sig")] = "purity"
        lnames[lnames.index("sigbulge")] = "mean_m"
        hcat.data.dtypes = tuple(lnames)
        
        #
        cal_purity(hcat)

        large_matched_gals = gal_hal_match(gcat, hcat, mgal_min=3e8)
        purity_crit = 0.995
        good_gals = gals_of_interest[gals_of_interest["host_purity"] > purity_crit]
        
        pickle.dump(good_gals, open("good_gals_{}.pickle".format(nout)))

