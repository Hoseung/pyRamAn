import numpy as np



def get_cell(allcell, region): 
    x=allcell["x"] 
    y=allcell["y"] 
    z=allcell["z"] 
 
    xc, yc, zc = region["centers"] 
    rr = region["radius"] 
 
    ind = np.where(np.square(x-xc) + 
                   np.square(y-yc) + 
                   np.square(z-zc) < rr)[0] 
 
    # advanced indexing = copy. 
    return allcell[ind] 
 
     
def get_cell_kd(kdtree, gal, rscale=15.0): 
    """ 
    Extract cells within rscale * Rreff and add to the galaxy.  
    """ 
    xc,yc,zc = gal.meta.xc, gal.meta.yc, gal.meta.zc 
    rgal = gal.meta.reff * rscale 
    __distance, index = kdtree.query_ball_point((xc,yc,zc), rgal) 
    gg.cell = kdtree.data[index] 
 

def do_work(sub_sample, nout,
            rscale=2.0):
    """ 
    Per process. 
    """ 
    import utils.sampling as smp 
    import galaxymodule as galm 
    from utils.cosmology import Timeconvert 
    from galaxymodule import mk_gal_params as mgp    
    from scipy.spatial import cKDTree
    from load.sim import Sim
    
    gen_vmap_sigmap_params = dict(npix_per_reff=5,
                                  rscale=3.0,
                                  n_pseudo=60,
                                  verbose=False,
                                  voronoi=None, #voronoi_dict
                                  weight="luminosity")


    cal_lambda_params = dict(npix_per_reff=5,
                             rscale=3.0,
                             method='ellip',
                             verbose=False,
                             iterate_mge = False,
                             save_result = True,
                             galaxy_plot_dir='./')
    mgp.HAGN["verbose"] = False
    Mcut = 1e10
 
    s = Sim(nout=nout)
    tc = timeconverter = Timeconvert(s.info)
    print("11")
    xrange = [min(sub_sample["x"] - sub_sample["r"] * rscale),
          max(sub_sample["x"] + sub_sample["r"] * rscale)]
    yrange = [min(sub_sample["y"] - sub_sample["r"] * rscale),
          max(sub_sample["y"] + sub_sample["r"] * rscale)]
    zrange = [min(sub_sample["z"] - sub_sample["r"] * rscale),
          max(sub_sample["z"] + sub_sample["r"] * rscale)]
    #print(xrange, yrange, zrange)
    region = smp.set_region(ranges=[xrange, yrange, zrange])


    s.set_ranges(region["ranges"])
    print("So far so good")
    print(s.ranges)
    return
    s.add_hydro(nvarh=5)
    
    
    kdtree = cKDTree(np.stack(s.hydro.cell["x"],
                              s.hydro.cell["y"],
                              s.hydro.cell["z"]))
    
    for gcat_this in sub_sample:
        gg = galm.rd_GM.Gal(nout=nout, 
                            catalog=gcat_this.copy(),
                            info=s.info)
        gg.debug=False
        make_gal.mk_gal(gg,**mk_gal_params)
        
        gg.star['time'] = tc.time2gyr(gg.star['time'],
                                        z_now = gg.info.zred)
        get_cell_kd(kdtree, gg)
        if save:
            pickle.dump(gg.cell, open("CELL_"+str(nout) + "_" + str(gg.meta.id) + ".pickle", "wb"))
        # convert to kpc
        gg.cell["x"] = (gg.cell["x"]-gg.center_code[0])*gg.info.boxtokpc
        gg.cell["y"] = (gg.cell["y"]-gg.center_code[1])*gg.info.boxtokpc
        gg.cell["z"] = (gg.cell["z"]-gg.center_code[2])*gg.info.boxtokpc
        gg.cell["dx"] *= gg.info.boxtokpc    
        # Do other calculations
        gg.mgas_tot = np.sum(gg.cell["var0"]*gg.cell["dx"]**3)
        make_gal.extract_cold_gas(gg, dr=5, rmax=200)

        # Now star and cell memberships are determined. 
    
        do_my_jobs(gg)


def domain_decompose_cat(gcat, nbins=5):
    """
        divide catalog into nbins**3 cubics and return...
        yield each chunk of cat.data.

        gcat is a partial catalog: ind != id -1
    """
    import utils.match as mtc
    ind_all = np.floor(gcat.data["x"]*nbins).astype(int) \
            + nbins * np.floor(gcat.data["y"]*nbins).astype(int) \
            + nbins**2*np.floor(gcat.data["z"]*nbins).astype(int)

    ind_sort = np.argsort(ind_all)

    sd = sorted_data = gcat.data[ind_sort]
    sorted_ind_all = ind_all[ind_sort]

    for i in range(nbins**3):
        #i_now = np.where(sorted_ind_all == i)[0]
        #sub_sample_ind = sd[i_now]["id"] -1
        #inds =  mtc.match_list_ind(gcat.data["id"], sd[i_now]["id"])
        yield gcat.data[mtc.match_list_ind(gcat.data["id"], sd[np.where(sorted_ind_all == i)[0]]["id"])]


def cat_only_relevant_gals(gcat, all_sample_ids, nout):
    import utils.match as mtc
    allgal_now = np.array(all_sample_ids[str(nout)])
    gcat.data = gcat.data[mtc.match_list_ind(gcat.data["id"], allgal_now)]


