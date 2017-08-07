import numpy as np

def ind_cell(allcell, region):
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


def ind_cell_kd(kdtree, gal, pboxsize, rscale=25.0):
    """
    Extract cells within rscale * Rreff and add to the galaxy.
    """
    xc,yc,zc = gal.meta.xc, gal.meta.yc, gal.meta.zc
    rgal = gal.meta.reff * rscale / (pboxsize*1e3) # kpc -> code unit
    #index = kdtree.query_ball_point((xc,yc,zc), rgal)
    xyzcen = (xc/pboxsize + 0.5,
              yc/pboxsize + 0.5,
              zc/pboxsize + 0.5)
    return kdtree.query_ball_point(xyzcen, rgal)


def get_cell(allcell, kdtree, gg, info):
    # Simple spherical cut.
    gg.cell=allcell[ind_cell_kd(kdtree, gg, info.pboxsize)]

    if len(gg.cell) > 1:
        #print(s.hydro.cell["x"].ptp())
        # convert to kpc
        #print("gg.center_code", gg.center_code)
        gg.cell["x"] = (gg.cell["x"]-gg.center_code[0])*info.boxtokpc
        gg.cell["y"] = (gg.cell["y"]-gg.center_code[1])*info.boxtokpc
        gg.cell["z"] = (gg.cell["z"]-gg.center_code[2])*info.boxtokpc
        gg.cell["dx"] *= info.boxtokpc


def add_output_containers(gg):
    gg.meta.sfr_results={"hist_dt":None, "hist_tmin":None, "hist_tmax":None, "hist":None, "sfr_dts":None, "sfrs":None, "area":None}
    #gg.meta.lambda_results={"lambda_results", }
    #gg.meta.mge_results={"mge_results":None}
    gg.meta.gas_results={"gas_results":None, "mgas_tot":None, "mgas_cold":None, "Ln_gas":None}
    gg.meta.vsig_results={"Vmax":None, "sigma":None, "V_sig":None}

def do_work(sub_sample, nout, i_subsample,
            rscale=1.5, save_cell=False):
    """
    Per process.
    """
    import utils.sampling as smp
    from galaxymodule import make_gal
    from utils.cosmology import Timeconvert
    from galaxymodule import mk_gal_params as mgp
    from galaxymodule import rotation_parameter, rd_GM
    from scipy.spatial import cKDTree
    from load.sim import Sim
    from galaxymodule import vmax_sig
    import pickle
    from galaxymodule.quick_mock import Simplemock
    from galaxymodule import gal_properties
    import time
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

    sfr_params = dict(hist_dt=0.1,
                      hist_tmin=0,
                      hist_tmax=None,
                      sfr_dts = [0.1, 0.5, 1.0])

    gas_params = dict(dr=5, rmax=200, density_ratio=1e-3)


    mgp.HAGN["verbose"] = False
    mgp.HAGN["mstar_min"] = 1e7
    out_dir = "./lambda_results/" + str(nout) +'/'

    # Common 1
    # simulation information
    s = Sim(nout=nout)
    xrange = [min(sub_sample["x"] - sub_sample["r"] * rscale),
          max(sub_sample["x"] + sub_sample["r"] * rscale)]
    yrange = [min(sub_sample["y"] - sub_sample["r"] * rscale),
          max(sub_sample["y"] + sub_sample["r"] * rscale)]
    zrange = [min(sub_sample["z"] - sub_sample["r"] * rscale),
          max(sub_sample["z"] + sub_sample["r"] * rscale)]
    #print(xrange, yrange, zrange)
    region = smp.set_region(ranges=[xrange, yrange, zrange])

    t0 = time.time()
    s.set_ranges(region["ranges"])
    #print(s.ranges)
    # Common2
    # Hydro cell data
    s.add_hydro(nvarh=5)
    t1 = time.time()
    print("Loading hydro took", t1 - t0)
    # Common 3
    # Cell KDTree
    kdtree = cKDTree(np.stack((s.hydro.cell["x"],
                               s.hydro.cell["y"],
                               s.hydro.cell["z"])).T)

    # Common 4
    # Stellar age converter
    tc = timeconverter = Timeconvert(s.info)

    # Common 5
    # Mock image generator
    MockSED = Simplemock()#repo=dfl.dir_repo+'sed/')

    result_sub_sample=[]
    print("{} galaxies in this sub ample".format(len(sub_sample)))
    for i, gcat_this in enumerate(sub_sample):
        gg = rd_GM.Gal(nout=nout,
                       catalog=gcat_this.copy(),
                       info=s.info)
        #print("s.info.pboxsize", s.info.pboxsize)
        gg.debug=False
        make_gal.mk_gal(gg,**mgp.HAGN)
        gg.meta.idx = gcat_this["idx"]
        #print(i,"-th, idx=", gg.meta.idx)
        if gg.meta.nstar < 60:
            gg.meta.mgas_cold = -1
            result_sub_sample.append(gg.meta)
            continue

        gg.meta.nout = nout
        gg.star['time'] = tc.time2gyr(gg.star['time'],
                                        z_now = gg.info.zred)

        # Add output dicts
        add_output_containers(gg)

        # gas properties
        get_cell(s.hydro.cell, kdtree, gg, s.info)
        #gg.cell = s.hydro.cell[ind_cell_kd(s.hydro.cell, kdtree, gg, s.info)]
        if len(gg.cell) > 1:
            if save_cell:
                pickle.dump(gg.cell, open("CELL_"+str(nout) + "_" + str(gg.meta.id) + ".pickle", "wb"))

            gal_properties.get_cold_cell(gg, s.info, **gas_params)
            gal_properties.get_gas_properties(gg, s.info)

        # Now star and cell memberships are determined.

        # Cell needed to calcluate gas attenuation.
        # r-band luminosity to be used as weights.
        gg.star.Flux_r= MockSED.get_flux(star=gg.star, filter_name='r')

       	gg.meta.mean_age = np.average(gg.star["time"], weights=gg.star["m"])
        gg.cal_norm_vec()

        gg.meta.rscale_lambda = gen_vmap_sigmap_params["rscale"]
        # Make pseudo particles. - memory usage!
        rotation_parameter.gen_vmap_sigmap(gg, **gen_vmap_sigmap_params)
        rotation_parameter.cal_lambda_r_eps(gg, **cal_lambda_params)

        # Calculate Vmax, Sig
        # get_vmax_sig uses vmap and sigmap from get_vmap_sigmap.
        # So if the maps were luminosity weighted, that also applies to the vmax and sig.
        vmax_sig.get_vmax_sig(gg, make_plot=False)

        # SFR
        # send the parameters to the begining.
        # Each function assumes various attributes from the gal object.
        # you may want to check them exist before/in each function.
        gal_properties.get_sfr_all(gg, **sfr_params)

        # Misc


        result_sub_sample.append(gg.meta)

    fout = out_dir + "result_sub_sample_{}_{}.pickle".format(nout, i_subsample)
    pickle.dump(result_sub_sample, open(fout, "wb"))
    print("Galaxy properties took", time.time() - t1 )


#########################################################

def domain_decompose_cat(gcat, nbins=5):
    """
        divide catalog into nbins**3 cubics.
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
        yield gcat.data[mtc.match_list_ind(gcat.data["id"], sd[np.where(sorted_ind_all == i)[0]]["id"])]


def cat_only_relevant_gals(gcat, all_sample_ids, nout):
    import utils.match as mtc
    allgal_now = np.array(all_sample_ids[str(nout)])
    gcat.data = gcat.data[mtc.match_list_ind(gcat.data["id"], allgal_now)]
