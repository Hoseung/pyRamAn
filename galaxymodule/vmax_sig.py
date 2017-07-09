import numpy as np

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
                 make_plot=True,
                 nreff=2.0):

    # get points near the major axis line
    n_pseudo = 30
    img_size = 40
    npix_per_reff = 5

    if not hasattr(gal.meta, "mge_result_list"):
        lambdas = cal_lambda_r_eps(gal, save_result=False, n_pseudo=n_pseudo,
                                   npix_per_reff=npix_per_reff)


    xinds = np.tile(np.arange(img_size) + 0.5, img_size)
    yinds = np.repeat(np.arange(img_size) + 0.5, img_size)

    ## Pixels on the major axis.

    # polynomial constants
    fit = gal.meta.mge_result_list[0]
    f_a = np.tan(fit["pa_rad"])
    f_b = -1
    f_c = fit["ycen"] - f_a*fit["xcen"]

    distance_to_line = np.abs(f_a*xinds + f_b*yinds + f_c)/np.sqrt(f_a**2 + f_b**2)
    i_ok = distance_to_line < 1

    if make_plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2)

        ax[0].imshow(gal.vmap, origin="lower")
        ax[0].scatter(xinds[i_ok], yinds[i_ok], color="g")

    v_good = gal.vmap.flatten()[i_ok]
    sig_good = gal.sigmap.flatten()[i_ok]
    d_good = np.sqrt(np.square(xinds[i_ok]-fit["xcen"]) + np.square(yinds[i_ok]-fit["ycen"]))


    d_sort = np.argsort(d_good)
    v_good = np.abs(v_good[d_sort])
    sig_good= sig_good[d_sort]
    d_good = d_good[d_sort]

    npoints = len(v_good)
    binsize = np.int(npoints/10)

    # smoothed arrays
    v_smooth=np.array([np.mean(v_good[i*binsize:(i+1)*binsize]) for i in range(10)])
    d_smooth=np.array([np.mean(d_good[i*binsize:(i+1)*binsize]) for i in range(10)])

    imax = np.argmax(v_smooth[d_smooth<nreff*npix_per_reff])

    # sigma
    star = gal.star
    star = gal.star[np.where((np.square(star['x']) +
                        np.square(star['y']) +
                        np.square(star['z'])) < gal.meta.reff**2)[0]]# in kpc unit

    sig = weighted_std(star["vz"], star['m'])

    vmax = v_smooth[imax]
    if make_plot:
        plt.scatter(d_good, v_good, color="r")
        plt.scatter(d_good, sig_good, color="b")
        plt.plot(d_smooth, v_smooth, color="r")
        plt.scatter(d_smooth[imax], v_smooth[imax], s=200, marker="^", color="r")
        plt.scatter(0, sig, s=200, marker="v", color='b')

        print("Vmax {:.2f}, sigma {:.2f}, v/sig {:.2f}".format(vmax, sig, vmax/sig))
        ax[1].set_aspect('auto')
        plt.savefig(str(gal.meta.id) + "_vel_curve" + str(n_pseudo) + ".png")

    gal.meta.vsig_results["Vmax"]= vmax
    gal.meta.vsig_results["sigma"] = sig
    gal.meta.vsig_results["V_sig"] = vmax/sig
