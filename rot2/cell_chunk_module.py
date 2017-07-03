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
    gg.cell=allcell[ind_cell_kd(kdtree, gg, info.pboxsiz)]

    if len(gg.cell) > 1:
        #print(s.hydro.cell["x"].ptp())
        # convert to kpc
        #print("gg.center_code", gg.center_code)
        gg.cell["x"] = (gg.cell["x"]-gg.center_code[0])*info.boxtokpc
        gg.cell["y"] = (gg.cell["y"]-gg.center_code[1])*info.boxtokpc
        gg.cell["z"] = (gg.cell["z"]-gg.center_code[2])*info.boxtokpc
        gg.cell["dx"] *= info.boxtokpc


class Containers():
    def __init__(self):
        from collections import namedtuple
        results_sfr=namedtuple("sfr_results",
                ["dt", "hist_tmin", "hist_tmax", "h", "sfr_dts", "sfrs", "area"])
        results_lambda=namedtuple("lambda_results", [])
        results_mge=namedtuple("mge_results", [])
        results_gas=namedtupe("gas_results",
             ["mgas_tot", "mgas_cold", "Ln_gas"])



def do_work(sub_sample, nout,
            rscale=2.0, save_cell=False):
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
                      sfr_dt = [0.1, 0.5, 1.0],
                      sfr_has_hist=True)

    mgp.HAGN["verbose"] = False
    mgp.HAGN["mstar_min"] = 1e7
    out_dir = "./lambda_results/"

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


    s.set_ranges(region["ranges"])
    print(s.ranges)
    # Common2
    # Hydro cell data
    s.add_hydro(nvarh=5)

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

    # Common 6
    # Output containers
    CT = Containers() # CT can be imported directly like astropy.cosmology.WMAP.

    result_sub_sample=[]

    for gcat_this in sub_sample:
        #print("cat ", sub_sample["x"].ptp())
        gg = rd_GM.Gal(nout=nout,
                            catalog=gcat_this.copy(),
                            info=s.info)
        #print("s.info.pboxsize", s.info.pboxsize)
        print("loading galaxy ")
        gg.debug=False
        make_gal.mk_gal(gg,**mgp.HAGN)
        if gg.meta.nstar < 60:
            gg.meta.mgas_cold = -1
            result_sub_sample.append(gg.meta)
            continue

        gg.star['time'] = tc.time2gyr(gg.star['time'],
                                        z_now = gg.info.zred)

        # gas properties
        get_cell_kd(s.hydro.cell, kdtree, gg, s.info)
        if len(gg.cell) > 1:
            if save_cell:
                pickle.dump(gg.cell, open("CELL_"+str(nout) + "_" + str(gg.meta.id) + ".pickle", "wb"))

            gas.get_gas_all(gg, CT, **gas_params) # CT = output container.

        # Now star and cell memberships are determined.

        # Cell needed to calcluate gas attenuation.
        # r-band luminosity to be used as weights.
        gg.star.Flux_r= MockSED.get_flux(star=gg.star, filter_name='r')

       	gg.meta.mean_age = np.average(gg.star["time"], weights=gg.star["m"])
        gg.cal_norm_vec()

        # Can I put 'rscale_lambda' into a namedtuple??
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
        gal_properties.get_sfr_all(gg, CT, **sfr_params)

        # Misc
        gg.meta.nout = nout
        #gg.meta.idx = this_sat["idx"]

        result_sub_sample.append(gg.meta)

    fout = out_dir + "result_sub_sample_" + str(nout) + "_from" + str(sub_sample[0]["id"]) + ".pickle"
    pickle.dump(result_sub_sample, open(fout, "wb"))


#########################################################

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

#def do_my_jobs(gg):

def ind_dense(cell, rmax = 200, dr = 5, rmin=1):
    """
        Measure radial profile and returns indices of cells inside r_min,
        where r_min is the local minima of radial MASS profile.
        -> should I use density profile instead?
    """
    from scipy.signal import argrelmin
    # radial profile.
    rr = np.sqrt(np.square(cell["x"])+\
                 np.square(cell["y"])+\
                 np.square(cell["z"]))

    i_sort = np.argsort(rr)
    r_sorted = rr[i_sort]
    mm = cell["dx"]**3 * cell["var0"]
    m_sorted = mm[i_sort]
    rmax = max([10, min([np.max(rr), rmax])])
    #print(rmax)

    #print("rmax now", rmax)
    # Note 1.
    # Depends on the cell resolution. How about 8 * dx_min?
    # Larger dx will count in small satellites,
    # while smaller dx will make the measurement sensitive to density fluctuations.
    nbins= int(rmax/dr)

    frequency, bins = np.histogram(r_sorted, bins = nbins, range=[0, rmax])
    bin_centers = bins[:-1] + 0.5 * dr # remove the rightmost boundary.

    m_radial = np.zeros(nbins)
    ibins = np.concatenate((np.zeros(1,dtype=int), np.cumsum(frequency)))

    for i in range(nbins):
        m_radial[i] = np.sum(m_sorted[ibins[i]:ibins[i+1]])
        # Check stellar surface density
        sig_at_r = m_radial[i]/(2 * np.pi * bin_centers[i] * dr)

    # Find local minimum
    # 1. If there is flat zeros, take the first zero.
    # If not, use scipy.argrelmin
    i_zero = np.argmax(m_radial==0)
    if i_zero > 0:
        ind_min = i_zero -1
    else:
        try:
            ind_min= min(argrelmin(m_radial)[0]) -1 # 1D array for 1D input.
        except:
            ind_min = -1
        #ind_min = ind_min[np.argmax(ind_min * dr > rmin)]* dr

    # Note 2.
    # If the minimum is farther than rmin=10kpc,
    # I assume that is correct.
    return rr < bin_centers[ind_min]


def rho_t_cut(cell, info, lose_cut=False):
    """
        Extract galactic cold gas following Torrey+12 criterion.
        Assume cells in the original (code) unit.
    """
    # Var0 in Msun h^2 kpc^-3 unit.
    kpc_in_cm = 3.08567758e21
    msun_in_g = 1.99e33
    gcc2this_unit = kpc_in_cm**3/msun_in_g

    return np.log10(cell["var4"]/cell["var0"]*info.unit_T2) < 6 + 0.25*np.log10((cell["var0"]*info.unit_d)*gcc2this_unit*1e-10)#
