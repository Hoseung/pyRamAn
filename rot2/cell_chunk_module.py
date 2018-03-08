import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


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
    rgal = min([100, max([30, gal.meta.reff * rscale])]) / (pboxsize*1e3) # kpc -> code unit
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
        gg.cell["var1"] = gg.cell["var1"]*info.kms
        gg.cell["var2"] = gg.cell["var2"]*info.kms
        gg.cell["var3"] = gg.cell["var3"]*info.kms
        gg.cell["dx"] *= info.boxtokpc


def add_output_containers(gg):
    gg.meta.sfr_results={"hist_dt":None, "hist_tmin":None, "hist_tmax":None, "hist":None, "sfr_dts":None, "sfrs":None, "area":None}
    #gg.meta.lambda_results={"lambda_results", }
    #gg.meta.mge_results={"mge_results":None}
    gg.meta.gas_results={"gas_results":None, "mgas_tot":None, "mgas_cold":None, "Ln_gas":None}
    gg.meta.vsig_results={"Vmax":None, "sigma":None, "V_sig":None}

def do_work(sub_sample, nout, i_subsample,
            out_base='/scratchb01/hoseung/newrun/',
            wdir = "./",
            rscale=1.5,
            do_cell=True,
            out_subdir = "lambda_results/",
            save_cell=True,
            do_plot=False):
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
    import galaxymodule.quick_mock as qmc
    from galaxymodule import gal_properties
    import time
    import os

    #voronoi_dict=dict(targetSN=20, quiet=True, plot=False)
    voronoi_dict = None
    gen_vmap_sigmap_params = dict(npix_per_reff=5,
                                  rscale=3.0,
                                  n_pseudo=60,
                                  verbose=False,
                                  voronoi=voronoi_dict,
                                  weight="luminosity")
    gen_vmap_sigmap_params2 = dict(npix_per_reff=5,
                                  rscale=3.0,
                                  n_pseudo=1,
                                  verbose=False,
                                  voronoi=voronoi_dict, #voronoi_dict
                                  weight="luminosity")

    cal_lambda_params = dict(npix_per_reff=5,
                             rscale=3.0,
                             method='ellip',
                             verbose=False,
                             voronoi=voronoi_dict,
                             iterate_mge = False,
                             save_result = True,
                             galaxy_plot_dir='./')


    sfr_params = dict(hist_dt=0.1,
                      hist_tmin=0,
                      hist_tmax=None,
                      sfr_dts = [0.1, 0.5, 1.0])

    gas_params = dict(dr=5, rmax=200, density_ratio=1e-3)


    mgp.HAGN["verbose"] = False
    mgp.HAGN["mstar_min"] = 1e8 * nout / 782
    out_dir = out_base+out_subdir+ str(nout) +'/'

    # Common 1
    # simulation information
    s = Sim(base=wdir, nout=nout)

    if do_cell:
        print("Do CELL")
        # If there are CELL_ files available, use them.
        for this_gal in sub_sample:
            if os.path.isfile(out_base+"CELL_{:05d}/CELL_{:d}_{:d}.pickle".format(nout,nout,this_gal["id"])):
                this_gal["level"] = 1234 # Overwrite level to mark CELL_ availability.
                #print("Has cell")
            else:
                #print("No cell")
                this_gal["level"] = -9 # Overwrite level to mark CELL_ availability.

        if sum(sub_sample["level"] < 0) > 0:
            print("Missing CELL files")
            print(sub_sample["id"][sub_sample["level"] < 0])
            # has something to load
            xrange = [min(sub_sample["x"][sub_sample["level"] < 0] -\
                          sub_sample["r"][sub_sample["level"] < 0] * rscale),
                  max(sub_sample["x"][sub_sample["level"] < 0] +\
                      sub_sample["r"][sub_sample["level"] < 0] * rscale)]
            yrange = [min(sub_sample["y"][sub_sample["level"] < 0] -\
                          sub_sample["r"][sub_sample["level"] < 0] * rscale),
                  max(sub_sample["y"][sub_sample["level"] < 0] +\
                      sub_sample["r"][sub_sample["level"] < 0] * rscale)]
            zrange = [min(sub_sample["z"][sub_sample["level"] < 0] -\
                          sub_sample["r"][sub_sample["level"] < 0] * rscale),
                  max(sub_sample["z"][sub_sample["level"] < 0] +\
                      sub_sample["r"][sub_sample["level"] < 0] * rscale)]

            region = smp.set_region(ranges=[xrange, yrange, zrange])

            s.set_ranges(region["ranges"])
        # Common2
        # Hydro cell data
            t0 = time.time()
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
    MockSED = qmc.Simplemock(info=s.info)#repo=dfl.dir_repo+'sed/')
    print(s.info)
    sdss_band = qmc.BandSDSS()

    result_sub_sample=[]
    print("{} galaxies in this sub sample".format(len(sub_sample)))
    for i, gcat_this in enumerate(sub_sample):
        gg = rd_GM.Gal(nout=nout,wdir=wdir,
                       catalog=gcat_this.copy(),
                       info=s.info, type_cell=do_cell)
        #print("s.info.pboxsize", s.info.pboxsize)
        gg.debug=False
        mgp.HAGN["method_cov"] = "close_member"
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

        if do_cell:
            # gas properties
            fn_cell=out_base+"CELL_{:05d}/CELL_{:d}_{:d}.pickle".format(nout,nout,gg.meta.id)
            if gcat_this["level"] ==1234:
                gg.cell = pickle.load(open(fn_cell, "rb"))
                # Temporary fix for early CELLS
                if gg.cell["var1"].ptp() < 10:
                    gg.cell["var1"] = gg.cell["var1"]*info.kms
                    gg.cell["var2"] = gg.cell["var2"]*info.kms
                    gg.cell["var3"] = gg.cell["var3"]*info.kms

                print("Load cell")
            else:
                get_cell(s.hydro.cell, kdtree, gg, s.info)
                print("Read from hydro")
            #gg.cell = s.hydro.cell[ind_cell_kd(s.hydro.cell, kdtree, gg, s.info)]
                if len(gg.cell) > 1:
                    if save_cell:
                        pickle.dump(gg.cell, open(fn_cell,"wb"))

                gal_properties.get_cold_cell(gg, s.info, **gas_params)
                gal_properties.get_gas_properties(gg, s.info)
        else:
            gg.cell = None
            # Now star and cell memberships are determined.

        # Now stars and cells are ready. Correct metallicity
        gg.star["metal"] *=4.08 - 0.21*s.info.zred - 0.11*s.info.zred**2


        # Cell needed to calcluate gas attenuation.
        # r-band luminosity to be used as weights.
        gg.star.Flux_u= MockSED.get_flux(star=gg.star, cell = gg.cell, info=s.info, filter_name='u')
        gg.meta.Mu = qmc.get_absolute_mag(gg.star.Flux_u, band=sdss_band, bandname="u")
        gg.star.Flux_g= MockSED.get_flux(star=gg.star, cell = gg.cell, info=s.info, filter_name='g')
        gg.meta.Mg = qmc.get_absolute_mag(gg.star.Flux_g, band=sdss_band, bandname="g")
        gg.star.Flux_r= MockSED.get_flux(star=gg.star, cell = gg.cell, info=s.info, filter_name='r')
        gg.meta.Mr = qmc.get_absolute_mag(gg.star.Flux_r, band=sdss_band, bandname="r")

       	gg.meta.mean_age = np.average(gg.star["time"], weights=gg.star["m"])

        # SFR
        # send the parameters to the begining.
        # Each function assumes various attributes from the gal object.
        # you may want to check them exist before/in each function.

        gg.cal_norm_vec()

        gg.meta.rscale_lambda = gen_vmap_sigmap_params["rscale"]
        # Make pseudo particles. - memory usage!
        rotation_parameter.gen_vmap_sigmap(gg, **gen_vmap_sigmap_params)
        rotation_parameter.cal_lambda_r_eps(gg, **cal_lambda_params)
        if do_plot:
            plt.rcParams["mathtext.fontset"] = "dejavusans"
            fig,axs = plt.subplots(2,2)
            im = axs[0,0].imshow(gg.mmap, norm=LogNorm())
            fig.colorbar(im,ax=axs[0,0])
            axs[0,0].set_title("Stellar density")
            im=axs[0,1].imshow(gg.vmap, cmap="RdYlBu")
            fig.colorbar(im,ax=axs[0,1], label=r"$kms^{-1}$")
            axs[0,1].set_title("Velocity")
            im=axs[1,0].imshow(gg.sigmap)
            fig.colorbar(im,ax=axs[1,0], label=r"$kms^{-1}$")
            axs[1,0].set_title("Velocity dispersion")
            axs[1,1].plot(np.arange(15)/5.,gg.meta.lambda_result_list[0])
            axs[1,1].set_title("Stellar density")

            axs[0,0].set_xticks([0,10,20,30,40])
            xticks = ["{:.1f}".format(x) for x in np.arange(-4,5,2)]
            axs[0,0].set_xticklabels(xticks)
            axs[0,0].set_ylabel("position [kpc]")
            axs[0,0].set_yticks(np.linspace(0,40,5))
            yticks = ["{:.1f}".format(y) for y in np.linspace(-gg.meta.rgal, gg.meta.rgal, num=5)]
            axs[0,0].set_yticklabels(yticks)

            axs[0,1].set_xticks([0,10,20,30,40])
            xticks = ["{:.1f}".format(x) for x in np.arange(-4,5,2)]
            axs[0,1].set_xticklabels(xticks)
            axs[0,1].set_ylabel("position [kpc]")
            axs[0,1].set_yticks(np.linspace(0,40,5))
            yticks = ["{:.1f}".format(y) for y in np.linspace(-gg.meta.rgal, gg.meta.rgal, num=5)]
            axs[0,1].set_yticklabels(yticks)

            axs[1,0].set_xticks([0,10,20,30,40])
            xticks = ["{:.1f}".format(x) for x in np.arange(-4,5,2)]
            axs[1,0].set_xticklabels(xticks)
            axs[1,0].set_ylabel("position [kpc]")
            axs[1,0].set_yticks(np.linspace(0,40,5))
            yticks = ["{:.1f}".format(y) for y in np.linspace(-gg.meta.rgal, gg.meta.rgal, num=5)]
            axs[1,0].set_yticklabels(yticks)

            axs[0,0].set_xlabel(r"R/R$_{eff}$")
            axs[0,1].set_xlabel(r"R/R$_{eff}$")
            axs[1,0].set_xlabel(r"R/R$_{eff}$")
            axs[1,1].set_xlabel(r"R/R$_{eff}$")
            #axs[0,0].set_ylabel(r"R/R$_{eff}$")
            axs[1,1].set_ylabel(r"$\lambda$")
            axs[1,1].set_ylim([0,1.0])
            axs[1,1].set_aspect(3.)
            axs[1,1].set_xticks([0,1,2,3])
            xticks = ["{:1d}".format(x) for x in [0,1,2,3]]
            axs[1,1].set_xticklabels(xticks)
            axs[1,1].text(10, 0.8, "{:.2e}M" + r"$\odot$".format(gg.meta.mstar))
            fig.suptitle("ID {},  z={:.1f}".format(gg.meta.id, gg.info.zred))
            plt.tight_layout()
            plt.savefig(out_dir + "lam_map_{}_{}.png".format(nout, gg.meta.id), dpi=200)
            plt.close()

        # Calculate Vmax, Sig
        # get_vmax_sig uses vmap and sigmap from get_vmap_sigmap.
        # So if the maps were luminosity weighted, that also applies to the vmax and sig.
        vmax_sig.get_vmax_sig(gg, gg.meta.vsig_results, make_plot=False, out_dir=out_dir)

        # B/T?
        # Could be done after
        if False:
            gal.reorient()
            gal.cal_b2t(ptype='star', disk_criterion="Scannapieco2009",bound_only=False)
            # once again with edge-on?
            gg.meta.vsig_results_edge={"Vmax":None, "sigma":None, "V_sig":None}
            rotation_parameter.gen_vmap_sigmap(gg, **gen_vmap_sigmap_params2)
            rotation_parameter.cal_lambda_r_eps(gg, **cal_lambda_params)
            vmax_sig.get_vmax_sig(gg,gg.meta.vsig_results_edge, make_plot=False, out_dir=out_dir)

        # Calculate SFR before making pseudo particls
        gal_properties.get_sfr_all(gg, **sfr_params)


        # Misc
        result_sub_sample.append(gg.meta)

    fout = out_dir + "result_sub_sample_{}_{}.pickle".format(nout, i_subsample)
    #print(out_dir, fout)
    pickle.dump(result_sub_sample, open(fout, "wb"))
    #print("Galaxy properties took {:.2f}s with {} particles".format(time.time() - t1, gg.meta.nstar))
    return gg


#########################################################

def domain_decompose_cat(gdata, nbins=5):
    """
        divide catalog into nbins**3 cubics.
        yield each chunk of cat.data.

        gcat is a partial catalog: ind != id -1
    """
    import utils.match as mtc
    ind_all = np.floor(gdata["x"]*nbins).astype(int) \
            + nbins * np.floor(gdata["y"]*nbins).astype(int) \
            + nbins**2*np.floor(gdata["z"]*nbins).astype(int)

    ind_sort = np.argsort(ind_all)

    sd = sorted_data = gdata[ind_sort]
    sorted_ind_all = ind_all[ind_sort]

    for i in range(nbins**3):
        yield gdata[mtc.match_list_ind(gdata["id"], sd[np.where(sorted_ind_all == i)[0]]["id"])]


def cat_only_relevant_gals(gdata, all_sample_ids, nout):
    import utils.match as mtc
    allgal_now = np.array(all_sample_ids[str(nout)])
    gdata = gdata[mtc.match_list_ind(gdata["id"], allgal_now)]
