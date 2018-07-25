import numpy as np
#from galaxymodule import make_gal
#from galaxymodule import mk_gal_params as mgp
#from galaxymodule import rd_GM
#from galaxymodule import rotation_parameter
#import utils
#from load import sim
#from scipy.spatial import cKDTree


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
        gg.cell["vx"] = gg.cell["vx"]*info.kms
        gg.cell["vy"] = gg.cell["vy"]*info.kms
        gg.cell["vz"] = gg.cell["vz"]*info.kms
        gg.cell["dx"] *= info.boxtokpc


def plot_rot_map(gg, npix_reff=5):
    rscale=gg.meta.rscale_lambda
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
    n_lam_points = len(gg.meta.lambda_result_list[0])
    axs[1,1].plot(np.arange(n_lam_points)/npix_reff,
                  gg.meta.lambda_result_list[0])
    axs[1,1].set_title("Stellar density")

    ticks = npix_reff * np.arange(0, 2*(rscale+1)+1, 2)
    axs[0,0].set_xticks(ticks)
    xticks_label = ["{:.1f}".format(x) for x in np.arange(-(rscale+1),rscale+2,2)]
    axs[0,0].set_xticklabels(xticks_label)
    axs[0,0].set_ylabel("position [kpc]")
    axs[0,0].set_yticks(ticks)
    yticks = ["{:.1f}".format(y) for y in np.linspace(-gg.meta.rgal, gg.meta.rgal, num=5)]
    axs[0,0].set_yticklabels(yticks)

    axs[0,1].set_xticks(ticks)
    axs[0,1].set_xticklabels(xticks_label)
    axs[0,1].set_ylabel("position [kpc]")
    axs[0,1].set_yticks(ticks)
    axs[0,1].set_yticklabels(yticks)

    axs[1,0].set_xticks(ticks)
    axs[1,0].set_xticklabels(xticks_label)
    axs[1,0].set_ylabel("position [kpc]")
    axs[1,0].set_yticks(ticks)
    axs[1,0].set_yticklabels(yticks)

    axs[0,0].set_xlabel(r"R/R$_{eff}$")
    axs[0,1].set_xlabel(r"R/R$_{eff}$")
    axs[1,0].set_xlabel(r"R/R$_{eff}$")
    axs[1,1].set_xlabel(r"R/R$_{eff}$")
    #axs[0,0].set_ylabel(r"R/R$_{eff}$")
    axs[1,1].set_ylabel(r"$\lambda$")
    axs[1,1].set_ylim([0,1.0])
    axs[1,1].set_aspect(rscale)
    axs[1,1].set_xticks([0,1,2,3])
    xticks = ["{:1d}".format(x) for x in [0,1,2,3]]
    axs[1,1].set_xticklabels(xticks)
    axs[1,1].text(10, 0.8, "{:.2e}M" + r"$\odot$".format(gg.meta.mstar))
    fig.suptitle("ID {},  z={:.1f}".format(gg.meta.id, gg.info.zred))
    plt.tight_layout()
    plt.savefig("lam_map_{}_{}_lum.png".format(nout, gg.meta.id), dpi=200)
    plt.close()


def add_output_containers(gg):
    gg.meta.sfr_results={"hist_dt":None,
                         "hist_tmin":None, "hist_tmax":None,
                         "hist":None, "sfr_dts":None, "sfrs":None, "area":None}
    #gg.meta.lambda_results={"lambda_results", }
    #gg.meta.mge_results={"mge_results":None}
    gg.meta.gas_results={"gas_results":None, "mgas_tot":None,
                         "mgas_cold":None, "Ln_gas":None}
    gg.meta.vsig_results={"Vmax":None, "sigma":None, "V_sig":None}
