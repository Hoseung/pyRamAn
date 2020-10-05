import numpy as np
from .rotation_parameter import cal_lambda_r_eps

def weighted_std(values, weights):
    import numpy as np
    import math
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise

    return math.sqrt(variance)

def get_vmax_sig(gal,
                 mge_results='',
                 make_plot=False,
                 nreff=2.0,
                 out_dir="./"):
    """
    Determine maximum rotation velocity and the projected velocity dispersion
    within 1Reff. Requires vmap and sigmap.

    parameters
    ----------
    nreff : 2.0
        maximum range in Reff where the Vmax is expected to occur.

    """

    npix_per_reff = gal.params.vmap_sigmap["npix_per_reff"]
    rscale = gal.params.vmap_sigmap["rscale"]

    # get points near the major axis
    results = getattr(gal.meta, mge_results)
    try:
        fit = results['mge_result_list'][0]
    except:
        print(f"Error... Can't find gal.meta.{mge_results}['mge_result_list']")
    

    img_size = gal.vmap.shape[0]
    xinds = np.tile(np.arange(img_size) + 0.5, img_size)
    yinds = np.repeat(np.arange(img_size) + 0.5, img_size)
    vmap = gal.vmap
    sigmap = gal.sigmap

    ## Pixels on the major axis.
    # polynomial constants
    fit = results['mge_result_list'][0]
    f_a = np.tan(fit["pa_rad"])
    f_b = -1
    f_c = fit["ycen"] - f_a*fit["xcen"]
    distance_to_line = np.abs(f_a*xinds + f_b*yinds + f_c)/np.sqrt(f_a**2 + f_b**2)
    i_ok = np.where(distance_to_line < 1)[0]

    if make_plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2)

        ax[0].imshow(vmap, origin="lower")
        ax[0].scatter(xinds[i_ok], yinds[i_ok], color="g")

    v_good = vmap.flatten()[i_ok]
    d_good = np.sqrt(np.square(xinds[i_ok]-fit["xcen"]) + np.square(yinds[i_ok]-fit["ycen"]))

    d_sort = np.argsort(d_good)
    v_good = np.abs(v_good[d_sort])
    d_good = d_good[d_sort]

    npoints = len(v_good)
    binsize = np.int(npoints/10)

    # smoothed arrays
    # necessary even with voronoi tesselation?
    v_smooth=np.array([np.mean(v_good[i*binsize:(i+1)*binsize]) for i in range(10)])
    d_smooth=np.array([np.mean(d_good[i*binsize:(i+1)*binsize]) for i in range(10)])

    imax = np.argmax(v_smooth[d_smooth<nreff*npix_per_reff])

    # sigma
    star = gal.star
    star = star[np.where((np.square(star['x']) +
                        np.square(star['y']) +
                        np.square(star['z'])) < gal.meta.reff**2)[0]]# in kpc unit

    sig = weighted_std(star["vz"], star['m'])

    vmax = v_smooth[imax]
    if make_plot:
        sig_good = sigmap.flatten()[i_ok]
        sig_good= sig_good[d_sort]

        plt.scatter(d_good, v_good, color="r")
        plt.scatter(d_good, sig_good, color="b")
        plt.plot(d_smooth, v_smooth, color="r")
        plt.scatter(d_smooth[imax], v_smooth[imax], s=200, marker="^", color="r")
        plt.scatter(0, sig, s=200, marker="v", color='b')

        print("Vmax {:.2f}, sigma {:.2f}, v/sig {:.2f}".format(vmax, sig, vmax/sig))
        ax[1].set_aspect('auto')
        plt.savefig(out_dir + "{}_{}_vel_curve.png".format(gal.nout, gal.meta.id))
        plt.close()

    #gal.meta.vsig_results = 
    return dict(Vmax=vmax, sigma=sig, V_sig=vmax/sig)


def get_vmax_sig_Cappellari2007(gal,
                                mge_results='',
                                voronoi=None,
                                make_plot=False,
                                rscale=1,
                                out_dir="./"):

    results = getattr(gal.meta, mge_results)
    try:
        fit = results['mge_result_list'][0]
    except:
        print(f"Error... Can't find gal.meta.{mge_results}['mge_result_list']")
        return False

    mge_params = results['mge_result_list'][0]
    #print('mge_params', mge_params)
    # mge_params should have ['xcen', 'ycen', 'cos', 'sin', 'sma', 'smi']
    npix_per_reff = gal.npix_per_reff

    if voronoi is not None:
        return _measure_vmax_sig_Cappellari2007(mge_params,
                       gal.mmap_v, gal.vmap_v, gal.sigmap_v,
                       gal.xNode, gal.yNode,
                       npix_per_reff,
                       rscale,
                       voronoi=True)
    else:
        xpos = None
        ypos = None
        return _measure_vmax_sig_Cappellari2007(mge_params,
                            gal.mmap.ravel(), gal.vmap.ravel(), gal.sigmap.ravel(),
                            xpos, ypos,
                            npix_per_reff,
                            rscale,
                            voronoi=False)


def _measure_vmax_sig_Cappellari2007(mge_par,
                    mmap, vmap, sigmap,
                    xNode, yNode,
                    npix_per_reff,
                    rscale,
                    voronoi=False,
                    do_plot=False,
                    verbose=False,
                    ):

    xcen=mge_par["xcen"]
    ycen=mge_par["ycen"]
    cos=mge_par["cos"]
    sin=mge_par["sin"]
    sma=mge_par["sma"]
    smi=mge_par["smi"]

    dd = np.sqrt(((xNode-xcen)*cos + (yNode-ycen)*sin)**2/sma**2 + \
                 ((yNode-ycen)*cos - (xNode-xcen)*sin)**2/smi**2) * \
                 npix_per_reff

    npix = round(npix_per_reff * rscale)
    #points = np.zeros(npix)

    if verbose: print("Reff = half light?1", sum(mmap[dd < 1.0])/ sum(mmap))
    dist1d = np.sqrt(np.square(xNode - xcen) + np.square(yNode - ycen))

    #for i in range(len(points)):
    #ind = np.where( (dd > i) & (dd < (i+1)))[0]
    ind = np.where(dd < (npix))[0]

    if len(ind) > 0:
        V = np.sum(mmap[ind] * vmap[ind]**2)
        if V > 0:
            ind2 = np.where(sigmap[ind] > 0)[0]
            Sig = np.sum(mmap[ind[ind2]] * sigmap[ind[ind2]]**2)

            return dict(Vmax=None, sigma=None, V_sig=V/Sig)
    return False



def _deprecated_get_vmax_sig_Cappellari2007(gal,
                 mge_results='',
                 voronoi=False,
                 make_plot=False,
                 out_dir="./"):
    """
    Determine maximum rotation velocity and the projected velocity dispersion
    within 1Reff following Cappellari2007. Requires vmap and sigmap.

    parameters
    ----------

    """

    npix_per_reff = gal.params.vmap_sigmap["npix_per_reff"]
    rscale = gal.params.vmap_sigmap["rscale"]
    # get points near the major axis
    results = getattr(gal.meta, mge_results)
    try:
        fit = results['mge_result_list'][0]
    except:
        print(f"Error... Can't find gal.meta.{mge_results}['mge_result_list']")
        return False

    if voronoi:
        vmap = gal.vmap_v
        sigmap = gal.sigmap_v
        mmap = gal.mmap_v
    else:
        vmap = gal.vmap
        sigmap = gal.sigmap
        mmap = gal.mmap

    img_size = vmap.shape[0]
    xinds = np.tile(np.arange(img_size) + 0.5, img_size)
    yinds = np.repeat(np.arange(img_size) + 0.5, img_size)
        
    ## Pixels within 1reff_elliptical
    # polynomial constants
    
    f_a = np.tan(fit["pa_rad"])
    f_b = -1
    f_c = fit["ycen"] - f_a*fit["xcen"]
    distance_eps = np.sqrt(np.square(xinds-fit["xcen"]) + np.square(yinds-fit["ycen"]))
    i_ok = np.where(distance_eps < (npix_per_reff))[0]

    v_good = vmap.flatten()[i_ok]
    sig_good = sigmap.flatten()[i_ok]
    flux_good = mmap.flatten()[i_ok]
    #dist_good = distance_eps[i_ok]

    v_sig = np.sum(flux_good*v_good**2)/np.sum(flux_good*sig_good**2)

    #gal.meta.vsig_results_Capp = 
    return dict(Vmax=None, sigma=None, V_sig=v_sig)
