import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


msun_in_g = 1.989e33
kpc_in_cm = 3.086e+21


def ind_cell_kd(kdtree, gal, pboxsize, rscale=25.0):
    """
    Extract cells within rscale * Rreff and add to the galaxy.
    """
    xc,yc,zc = gal.meta.xc, gal.meta.yc, gal.meta.zc
    rgal = min([100, max([30, gal.meta.reff * rscale])]) / (pboxsize*1e3)
    # kpc -> code unit
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

def plot_2d_simple_maps(gg):
    fig, axs = plt.subplots(2,2)
    axs = axs.ravel()
    fig.set_size_inches(12,8)
    axs[0].hist(gg.star["metal"],
                weights= gg.star["m"],
                histtype="step",
                 label="star")
    axs[0].legend()
    axs[0].set_xlabel("Z")

    axs[1].hist2d(gg.star["z"], gg.star["x"],
                  range=[[-25,25],[-25,25]], bins=100, norm=LogNorm())
    axs[2].hist2d(gg.star["y"], gg.star["z"],
                  range=[[-25,25],[-25,25]], bins=100, norm=LogNorm())
    axs[3].hist2d(gg.star["x"], gg.star["y"],
                  range=[[-25,25],[-25,25]], bins=100, norm=LogNorm())
    axs[1].set_aspect("equal")
    axs[2].set_aspect("equal")
    axs[3].set_aspect("equal")


def plot_2d_simple_gas_maps(gg):
    cell_to_msun = gg.info.unit_d/msun_in_g * kpc_in_cm**3
    cell_mass = gg.cell["rho"]*gg.cell["dx"]**3 * cell_to_msun

    # 2d maps
    fig, axs = plt.subplots(2,2)
    axs = axs.ravel()
    rgal = gg.meta.reff * 5
    axs[0].hist2d(gg.star["x"], gg.star["y"], weights=gg.star["m"],
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm())
    axs[1].hist2d(gg.dm["x"], gg.dm["y"], weights=gg.dm["m"],
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm())
    axs[2].hist2d(gg.cell["x"], gg.cell["y"], weights=cell_mass,
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm(), vmin=1e4)
    for ax in axs:
        ax.set_aspect("equal")
    plt.savefig("{}_{}_gas.png".format(gg.nout, gg.meta.id), dpi=200)


def plot_radial(gg, nbins=50, rmax=3):
    gg.info.unit_d/msun_in_g * kpc_in_cm**3
    cell_to_msun = gg.info.unit_d/msun_in_g * kpc_in_cm**3

    fig, ax = plt.subplots()
    bins = np.linspace(0,3,50)**3

    dist_star = np.sqrt(np.sum(np.square(gg.star["pos"]), axis=1))
    isort_dist_star = np.argsort(dist_star)

    # star
    h_st, bin_st = np.histogram(dist_star[isort_dist_star], bins=bins,
                         weights=gg.star["m"][isort_dist_star])
    ax.plot(bins[1:]**1/3, h_st, label="star")

    if hasattr(gg, "dm") and gg.dm is not None:
        dist_dm = np.sqrt(np.sum(np.square(gg.dm["pos"]), axis=1))
        isort_dist_dm = np.argsort(dist_dm)

        h_dm, bin_dm = np.histogram(dist_dm[isort_dist_dm],  bins=bins,
                                    weights=gg.dm["m"][isort_dist_dm])
        ax.plot(bins[1:]**1/3, h_dm, label="dm")

    if hasattr(gg, "cell") and gg.cell is not None:
        dist_cell = np.sqrt(np.sum(np.square(gg.cell["pos"]), axis=1))
        isort_dist_cell = np.argsort(dist_cell)

        h_cell, bin_cell = np.histogram(dist_cell[isort_dist_cell],  bins=bins,
                 weights=gg.cell["rho"][isort_dist_cell]*
                         gg.cell["dx"][isort_dist_cell]**3 * cell_to_msun)
        ax.plot(bins[1:]**1/3, h_cell, label="gas")

    # all
    all_comp = h_st + h_dm + h_cell
    ax.plot(bins[1:]**1/3, h_st + h_dm + h_cell,
            label="all")

    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0.1,3**3])
    ax.set_ylim([1e7,1e11])
    ax.set_xlabel("kpc")
    ax.set_ylabel(r"$M_{\odot}$")
    plt.savefig("{}_{}_comp_profile.png".format(gg.nout, gg.meta.id), dpi=200)


def plot_bulge_disk(gg, bulge, disk, nbins=200, band="flux_u"):
    fig, axs = plt.subplots(3,2)
    fig.set_size_inches(6,9)
    axs = axs.ravel()
    axs[0].hist2d(bulge["x"], bulge["y"],
                  weights=bulge[band], bins=nbins, norm=LogNorm())
    axs[2].hist2d(bulge["y"], bulge["z"],
                  weights=bulge[band], bins=nbins, norm=LogNorm())
    axs[4].hist2d(bulge["z"], bulge["x"],
                  weights=bulge[band], bins=nbins, norm=LogNorm())
    axs[1].hist2d(disk["x"], disk["y"],
                  weights=disk[band], bins=nbins, norm=LogNorm())
    axs[3].hist2d(disk["y"], disk["z"],
                  weights=disk[band], bins=nbins, norm=LogNorm())
    axs[5].hist2d(disk["z"], disk["x"],
                  weights=disk[band], bins=nbins, norm=LogNorm())
    axs[0].set_title("Bulge  (e < 0.5)")
    axs[1].set_title("disk  (e > 0.8)")
    for ax in axs:
        ax.set_aspect("equal")

    fig.suptitle("B/T (0.5) = {:.2f}".format(np.sum(bulge["m"])/np.sum(gg.star["m"])))
    plt.savefig("{}_{}_xyz_{}-band.png".format(gg.nout, gg.meta.id, band),
                dpi=200)


def plot_rot_map(gg):
    """
    Make a rotation parameter with 4 panels.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    npix_reff = gg.params.vmap_sigmap["npix_per_reff"]
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
    plt.savefig("lam_map_{}_{}_lum.png".format(gg.nout, gg.meta.id), dpi=200)
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
