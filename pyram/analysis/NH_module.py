import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

msun_in_g = 1.989e33
kpc_in_cm = 3.086e+21

def load_components(gg, s, idlist=None,
                    load_dm=True,
                    load_cell=True,
                    load_raw=False,
                    save_cell=False,
                    verbose=False,
                    gas_radius=None):
    """
    Parameters
    ----------
    gg : Galaxy instance
    s  : Sim instance
    idlist : member DM particle ID list
        If present, member DM particles are determined accordingly.
    load_dm :
    load_cell :
    load_raw : False
        If True, read hydro cell from the raw data.
        If False, load pickled cell.
    save_cell : False
        If True, pickle current cell data.
        If load_raw is True, then save_cell is overidden to True.
    gas_radius : None, [in code unit]
        Specify the radisu of gas cell.
        By default, all gas inside 1Rvir is considered,
        but sometimes it is interesting to consider gas that are expelled out to ~3Rvir.
    """
    import pickle
    from ..utils import match as mtc
    from ..utils.sampling import Region
    from scipy.spatial import cKDTree
    nout = s.nout

    if load_dm or load_cell:
        reg = Region()
        reg.region_from_halo(gg.hcat)
    else:
        return

    if load_dm:
        s.set_ranges(reg.ranges)
        s.add_part(ptypes=["dm id pos vel mass"])
        ind = mtc.match_list_ind(s.part.dm["id"], idlist, allow_swap=False)
        gg.dm = s.part.dm[ind]
        gg.dm["pos"] -= gg.center_code
        gg.dm["pos"] *= gg.info.boxtokpc
        gg.dm["vel"] *= gg.info.kms
        gg.dm["m"] *= gg.info.msun

    if load_cell:
        # gas properties
        fn_cell=s.base+"/GalaxyMaker/CELL_"+\
                "{:05d}/CELL_{:d}_{:d}.pickle".format(nout,nout,gg.meta.id)

        if not load_raw:
            try:
                gg.cell = pickle.load(open(fn_cell, "rb"))
            except:
                pass
        else:
            reg.radius = max([reg.radius, 1e-5])
            if load_raw: save_cell = True

            if gas_radius is not None:
                gas_reg = Region()
                gas_reg.region_from_halo(gg.hcat)
                gas_reg.radius = gas_radius
                # Use region and update radius.
            else:
                # Use the same region as stars
                gas_reg = reg

            gg.gas_region = gas_reg # Note that gg.region is in kpc unit.
            s.set_ranges(gg.gas_region.ranges)
            s.add_hydro(verbose=False, load=True, pure=False)
            kdtree = cKDTree(np.stack((s.hydro.cell["x"],
                                       s.hydro.cell["y"],
                                       s.hydro.cell["z"])).T)
            get_cell(s.hydro.cell, kdtree, gg, gas_radius=gas_radius)
            #print("Retrieving cells from sim.hydro")
            if len(gg.cell) > 1 and save_cell:
                pickle.dump(gg.cell, open(fn_cell,"wb"))


def ind_cell_kd(kdtree, gal, pboxsize, gas_radius="25Reff"):
    """
    Extract cells within gas_radius and add to the galaxy.

    Todo:
    parsing part should appear earlier in load_components.
    """
    xc,yc,zc = gal.meta.xc, gal.meta.yc, gal.meta.zc
    if gas_radius is not None:
        if isinstance(gas_radius, str):
            if "Reff" in gas_radius:
                rscale_reff = float(gas_radius.split("Reff")[0])
                rgas = rscale_reff * gal.meta.reff
            elif "kpc" in gas_radius:
                rgas = float(gas_radius.split("kpc")[0])/pboxsize*1e3
            elif "Mpc" in gas_radius:
                rgas = float(gas_radius.split("Mpc")[0])/pboxsize
        elif isinstance(gas_radius, (int, float)):
            rgas = gas_radius
    return kdtree.query_ball_point((xc,yc,zc), rgas)


def get_cell(allcell, kdtree, gg, gas_radius="25Reff"):
    """
    It is non-trivial
    """
    # Simple spherical cut.
    gg.cell=allcell[ind_cell_kd(kdtree,
                                gg,
                                gg.info.pboxsize,
                                gas_radius=gas_radius)]

    if len(gg.cell) > 1:
        gg.cell["x"] = (gg.cell["x"]-gg.center_code[0])*gg.info.boxtokpc
        gg.cell["y"] = (gg.cell["y"]-gg.center_code[1])*gg.info.boxtokpc
        gg.cell["z"] = (gg.cell["z"]-gg.center_code[2])*gg.info.boxtokpc
        gg.cell["vx"] = gg.cell["vx"]*gg.info.kms
        gg.cell["vy"] = gg.cell["vy"]*gg.info.kms
        gg.cell["vz"] = gg.cell["vz"]*gg.info.kms
        gg.cell["dx"] *= gg.info.boxtokpc

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
    plt.savefig(f"{gg.nout}_{gg.meta.id}_simple.png", dpi=200)
    plt.close()

def plot_2d_simple_gas_maps(gg):
    cell_to_msun = gg.info.unit_d/msun_in_g * kpc_in_cm**3
    cell_mass = gg.cell["rho"]*gg.cell["dx"]**3 * cell_to_msun

    # 2d maps
    fig, axs = plt.subplots(3,3)
    axs = axs.ravel()
    rgal = gg.meta.reff * 5
    axs[0].hist2d(gg.star["x"], gg.star["y"], weights=gg.star["m"],
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm())
    axs[1].hist2d(gg.star["x"], gg.star["y"], weights=gg.star["m"],
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm())
    axs[2].hist2d(gg.star["x"], gg.star["y"], weights=gg.star["m"],
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm())
    axs[3].hist2d(gg.dm["x"], gg.dm["y"], weights=gg.dm["m"],
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm())
    axs[4].hist2d(gg.dm["x"], gg.dm["y"], weights=gg.dm["m"],
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm())
    axs[5].hist2d(gg.dm["x"], gg.dm["y"], weights=gg.dm["m"],
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm())
    axs[6].hist2d(gg.cell["x"], gg.cell["y"], weights=cell_mass,
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm(), vmin=1e4)
    axs[7].hist2d(gg.cell["x"], gg.cell["y"], weights=cell_mass,
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm(), vmin=1e4)
    axs[8].hist2d(gg.cell["x"], gg.cell["y"], weights=cell_mass,
                  bins=100, range=[[-rgal,rgal]]*2, norm=LogNorm(), vmin=1e4)
    for ax in axs:
        ax.set_aspect("equal")
    plt.savefig(f"{gg.nout}_{gg.meta.id}_gas.png", dpi=200)
    plt.close()

def cal_radial(data, nbins=50, rmax=25,
                surface_density=False,
                cell_to_msun=None):
    """
    return radial profile of the given species.
    Assumes to be centered at [0,0,0].

    parameters
    ----------
    nbins :
    rmax : 25kpc.
        Maximum extent of the galaxy in kpc unit.
    surface_density : False
        return 2D projected value rather than 3D density.
    cell_to_msun : None
        If RAMSES cell, a conversion factor is needed.
    """
    bins = np.linspace(0,rmax,nbins)

    dist = np.sqrt(np.sum(np.square(data["pos"]), axis=1))
    isort_dist = np.argsort(dist)

    # star
    try:
        h, hbin = np.histogram(dist[isort_dist], bins=bins,
                             weights=data["m"][isort_dist])
    except:
        assert cell_to_msun is not None, "Need cell_to_msun!"
        h, hbin = np.histogram(dist[isort_dist],  bins=bins,
                 weights=data["rho"][isort_dist]*
                         data["dx"][isort_dist]**3 * cell_to_msun)

    bin_centers = 0.5*(hbin[:-1] + hbin[1:])
    bin_widths  = hbin[1:] - hbin[:-1]

    if surface_density:
        two_pi_r_dr = 2 * np.pi * bin_centers * bin_widths
        h /= two_pi_r_dr
    else:
        #print(bin_centers)
        #print(bin_widths)
        four_third_pi_r_r_dr = 4/3 * np.pi * (hbin[1:]**2 - hbin[:-1]**2)
        h /= four_third_pi_r_r_dr
    return hbin, h

def plot_radial(gg, nbins=50, rmax=25, xlog=False, ylog=True,
                surface_density=False):
    """
    parameters
    ----------
    surface_density : boolean
        If True, the radial profile is divided by r^2.
    """
    fig, ax = plt.subplots()
    bins, h_st = cal_radial(gg.star,
                             nbins=50,
                             rmax=25,
                             surface_density=surface_density)
    ax.plot(bins[1:], h_st, label="star")
    h_all = h_st

    if hasattr(gg, "dm") and gg.dm is not None:
        bins, h_dm = cal_radial(gg.dm,
                                 nbins=50,
                                 rmax=25,
                                 surface_density=surface_density)
        ax.plot(bins[1:], h_dm, label="dm")
        h_all += h_dm
    else:
        h_dm = []

    if hasattr(gg, "cell") and gg.cell is not None:
        cell_to_msun = gg.info.unit_d/msun_in_g * kpc_in_cm**3
        bins, h_cell = cal_radial(gg.cell,
                                nbins=50,
                                rmax=25,
                                surface_density=surface_density,
                                cell_to_msun=cell_to_msun)
        ax.plot(bins[1:], h_cell, label="gas")
        h_all += h_cell
    else:
        h_cell = []

    # all
    ax.plot(bins[1:], h_all,
            label="all")

    ax.legend()
    if xlog:
        ax.set_xscale("log")
        ax.set_xlim([0.1,rmax])
    else:
        ax.set_xlim([0.01,rmax])
    if ylog:
        ax.set_yscale("log")
        if surface_density:
            ax.set_ylim([5e5,5e10])
        else:
            ax.set_ylim([5e5,5e10])
    else:
        ax.set_ylim([1e8,1e10])
    if surface_density:
        ylabel = r"$\Sigma \, [M_{\odot}/kpc^2]$"
    else:
        ylabel = r"$M_{\odot} kpc^{-3}$"
    ax.set_xlabel("kpc")
    ax.set_ylabel(ylabel)
    if surface_density:
        fig_type = "surface density"
    else:
        fig_type = "density"
    plt.savefig("{}_{}_comp_{}_profile.png".format(gg.nout, gg.meta.id, fig_type),
                dpi=200)
    plt.close()

    return dict(bins=bins,
                star=h_st,
                dm=h_dm,
                cell=h_cell)


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
    plt.close()

def plot_rot_map(gg):
    """
    Make a rotation parameter plot with 4 panels.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    npix_reff = gg.params.vmap_sigmap["npix_per_reff"]
    rscale=gg.params.vmap_sigmap["rscale"]
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
    axs[1,1].set_title(r"$\lambda_{<r}$")

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


def get_j_profile(gg, binwidth="reff"):
    """
    Caculate the specific angular momentum profile of each component up to
    Rgal (which is determined by measuring the surface brightness of stars)
    in bins with a width of Reff.
    """
    from collections import namedtuple

    if binwidth=="reff":
        binwidth = gg.meta.reff
    mod = gg.meta.rgal//binwidth
    bins = np.arange(mod+2)*binwidth

    _, j_spec_s_pro, m_prof_s = cal_j_profile(gg.star, bins=bins)
    _, j_spec_d_pro, m_prof_d = cal_j_profile(gg.star[gg.star["ellip"] > gg.params.decompose["e_disk"][0]], bins=bins)
    _, j_spec_b_pro, m_prof_b = cal_j_profile(gg.star[gg.star["ellip"] < gg.params.decompose["e_bulge"][1]], bins=bins)
    _, j_spec_g_pro, m_prof_g = cal_j_profile(gg.cell, info=gg.info, bins=bins)
    _, j_spec_dm_pro, m_prof_dm = cal_j_profile(gg.dm, bins=bins)

    Prof = namedtuple("profile", ["bins",
                                  "M_s", "M_d", "M_b", "M_dm", "M_g",
                                 "j_s", "j_d", "j_b", "j_dm", "j_g"])
    #Massj = namedtuple("")
    gg.profile = Prof(bins, m_prof_s, m_prof_d, m_prof_b, m_prof_dm, m_prof_g,
                   j_spec_s_pro, j_spec_d_pro, j_spec_b_pro, j_spec_dm_pro, j_spec_g_pro)

    gg.meta.MJ = {"j_s":j_spec_s_pro[-1],
                  "j_profile_s":j_spec_s_pro,
                  "j_b" :j_spec_b_pro[-1],
                  "j_profile_b":j_spec_b_pro,
                  "j_g" :j_spec_g_pro[-1],
                  "j_profile_g":j_spec_g_pro,
                  "j_d" :j_spec_d_pro[-1],
                  "j_profile_d":j_spec_d_pro,
                  "j_dm":j_spec_dm_pro[-1],
                  "j_profile_dm":j_spec_dm_pro}


def cal_j_profile(species, info=None, bins=10):
    """
    Calculate specific angular momentum profiles in bins.
    Returns cumulative j profile (for r' < r), and mass profile in each shell.
    """
    from scipy.stats import binned_statistic

    dist = np.sqrt(np.einsum("...i,...i",species["pos"],species["pos"]))
    RV = np.cross(species["pos"], species["vel"])

    try:
        cell_to_msun = info.unit_d/info.msun_in_g * info.kpc_in_cm**3
        mm = species["rho"] * species["dx"]**3 * cell_to_msun
    except:
        mm = species["m"]

    jj_tot = mm * np.sqrt(np.einsum("...i,...i",RV,RV))

    RVx_profile = binned_statistic(dist, RV[:,0]*mm, 'sum', bins=bins)[0]
    RVy_profile = binned_statistic(dist, RV[:,1]*mm, 'sum', bins=bins)[0]
    RVz_profile = binned_statistic(dist, RV[:,2]*mm, 'sum', bins=bins)[0]
    mm_profile, bin_edges, _ = binned_statistic(dist, mm, 'sum', bins=bins)

    jj_cum_profile = np.sqrt(np.square(np.cumsum(RVx_profile)) +
                             np.square(np.cumsum(RVy_profile)) +
                             np.square(np.cumsum(RVz_profile))) / np.cumsum(mm_profile)

    return bin_edges, jj_cum_profile, mm_profile
