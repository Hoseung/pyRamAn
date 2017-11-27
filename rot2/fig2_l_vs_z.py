## Fig2
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP7, z_at_value
from utils import match
from scipy.stats import gaussian_kde


def age2zred(lts):
    import astropy.units as u
    zreds = [z_at_value(WMAP7.age, ll * u.Gyr) for ll in lts]
    return zreds

def zreds2nouts(nouts, nnza):
    pass

def modify_ticks1(ax, nnza, nouts, nbins = 40):
    """
    Parameters
    ----------
    nbins:
        # of bins along the y-axis.
    """

    #aexps = nnza.a2b(nouts, "nout", "aexp")
    zreds = nnza.a2b(nouts, "nout", "zred")
    nnouts = len(nouts)

    # For a given list of nouts,
    # calculate a nice-looking set of zreds AND lookback times
    zz_target = [0.0, 0.5, 1.0, 2.0]
    #zz_target = [0.0, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0]
    # If 3 of 4 are matched, the 3 must be 0, 0.5, 1.0
    # So, later epohc should come first and be matched first.
    # If zz_target = [2.0, 1.0, 0.5, 0.0] and the nout spans
    # 300 ~ 787, 2.0, 1.0, and 0.5 appears as the ticks at wrong places.

    x_tick_pos = match.match_list_ind(nouts, np.array([782, 535, 343, 197]))# - nouts[-1]
    ax.set_xlim([0, nnouts])

    ax.set_xticks(x_tick_pos)#[::-1]
    ax.set_xticklabels(labels = ["{:0.1f}".format(z) for z in zz_target])


def modify_ticks2(ax2, nnza, nouts, nbins = 40):
    """
    Parameters
    ----------
    nbins:
        # of bins along the y-axis.
    """
    #aexps = nnza.a2b(nouts, "nout", "aexp")
    zreds = nnza.a2b(nouts, "nout", "zred")
    nnouts = len(nouts)

    # For a given list of nouts,
    # calculate a nice-looking set of zreds AND lookback times
    #zz_target = [2.0, 1.0, 0.5, 0.0]
    #x_tick_pos = match.match_list_ind(nouts[::-1], np.array([197, 343, 535, 782]))
    #nnouts = len(nouts)

    # For a given list of nouts,
    # calculate a nice-looking set of zreds AND lookback times

    # searchsorted requires arrays be in ascending order.

    #  at z = 3, lbt = 11.5243, age of universe = 2.18844
    lbt_targets=[0,5,8,12]
    nnza.a2b(lbt_targets, "age", "zred")
    # INCOMPLETE
    z_targets_u_age = age2zred(lbt_targets)
    print(z_targets_u_age)
    x_tick_pos = len(nouts) - np.searchsorted(zreds, z_targets_u_age)#[::-1]# + nout_ini# + nout_min
    print(x_tick_pos)
    #u_age_pos = np.searchsorted(zreds, )

    ax2.set_xticks(x_tick_pos)
    ax2.set_xticklabels(labels = ["{:.0f}".format(l) for l in u_age_targets])
    ax2.set_xlabel("Age of the universe (Gyr)")

def density_map(x, y, bandwidth=3., xbins=100j, ybins=80j, **kwargs):
    from sklearn.neighbors import KernelDensity
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[x.min():x.max():xbins,
                      y.min():y.max():ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)



def density_map_slow(x, y, sort=True):
    from scipy.stats import gaussian_kde
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    z /= max(z)

    idx = z.argsort()
    xx, yy = x[idx], y[idx]
    z = z[idx]

    #im = ax.scatter(xx, yy, c=z, s=50, edgecolor='')
    return xx,yy,z


def plot_lambda_evol(alldata, nouts,
                     nnza, ax=None,
                     data_type="fine",
                     density = "kernel_column",
                     add_errorbar = True,
                     fname=None,
                     cmap ="jet",
                     nhexbin=45,
                     nout_max=782):
    """
    im= plot_lambda_evol(alldata, nouts,
                 nnza_cell,
                 density = "kernel",
                 fname="./RUN2/figs/Lambda_evol.png",
                 cmap ="jet")

    nout_max = 782 : extrapolated values (783 ~ 787) are unreliable.
    """

    # compile data
    nnouts = len(nouts)
    ngals_tot = len(alldata)
    nout_min = min(nouts)

    #nstep_max = nnza.nnza["nstep"].max()
    #nstep_min = nstep_max - nnouts
    # Starting nouts of galaxies are different.
    if data_type == "fine":
        lambda_evol_all = np.zeros([ngals_tot, nnouts])
        for igal, gal in enumerate(alldata):
            if len(gal.finedata) == 0:
                continue
            if True:
                #nstep = gal.finedata["nstep"]
                ind = np.where((gal.finedata["nout"] > nout_min)*\
                               (gal.finedata["nout"] < nout_max))[0]
                lambda_evol_all[igal,ind] = gal.finedata["lambda_r"][ind]
            else:
                pass

        # Cut out unreliable data
        ind_good_nout = np.where(nouts <= nout_max)[0]
        lambda_evol_all = lambda_evol_all[:,ind_good_nout]
        nouts = nouts[ind_good_nout]
        nnouts = len(nouts)
        #print(lambda_evol_all)
    elif data_type == "coarse":
        lambda_evol_all = np.zeros([ngals_tot, nnouts])
        for igal, gal in enumerate(alldata):
            lambda_evol_all[igal][:] = gal.main_data["lambda_r"][:nnouts]

    elif data_type == "array":
        lambda_evol_all = alldata["lambda_r"]

    # change ticks
    #aexps = nnza.a2b(nouts, "nout", "aexp")
    zreds = nnza.a2b(nouts, "nout", "zred")

    if ax is None:
        fig, ax = plt.subplots()
        fig.set_size_inches(6,4)
        #fname=None
    #ax2 = ax.twiny()

    modify_ticks1(ax, nnza, nouts)
    #modify_ticks2(ax2, nnza, nouts)

    if density == "heat":
        nbins=100
        den_map = np.zeros((nbins, nnouts))
        for i in range(nnouts):
            #print(i)
            den_map[:,i], ypoints = np.histogram(lambda_evol_all[:,i], bins=nbins, range=[0,0.9])
            den_map[:,i] /= den_map[:,i].max()
        im = ax.imshow(den_map, origin="lower",
                       cmap=cmap,
                       interpolation="gaussian")
        ax.set_aspect('auto')

    elif density == "hexbin":
        xx = np.tile(np.arange(nnouts), ngals_tot)
        all_data = lambda_evol_all.ravel()
        ind_ok = all_data > 0.01
        #xx,yy,z = density_map(xx[ind_ok], all_data[ind_ok])
        im = ax.hexbin(xx[ind_ok], 10 * all_data[ind_ok],
                       gridsize=nhexbin,
                       bins=None,
                       cmap=cmap)

        ax.set_ylim([-0.5, 9])
        #fig = plt.gcf()
        #cb = fig.colorbar(im, ax=ax)
        #cb.set_label('log10(N)')
        y_tick_pos=np.arange(10)
        ax.set_yticks(y_tick_pos)#[::-1]
        ax.set_yticklabels(labels = ["{:0.1f}".format(0.1*y) for y in y_tick_pos])

    elif density == "kernel":
        xx = np.tile(np.arange(nnouts), ngals_tot)
        all_data = lambda_evol_all.ravel()
        ind_ok = all_data > 0.01
        xx,yy,zz = density_map(xx[ind_ok], all_data[ind_ok])
        #im = ax.scatter(xx, yy, c=z, s=50, edgecolor='', cmap=cmap)
        im = ax.pcolormesh(xx, yy, zz, cmap=cmap)
        lambda_range=[0.01, 0.8]
        yticks_ok=[0.0, 0.2, 0.4, 0.6, 0.8]
        ax.set_ylim([-0.05, 0.9])
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([str(yy) for yy in yticks_ok])
    elif density == "kernel_column":
        from matplotlib.patches import Circle
        #xx = np.tile(np.arange(nnouts), ngals_tot)
        patches = []
        all_xxs=[]
        all_yys=[]
        all_zzs=[]
        for inout, nout in enumerate(nouts):
            xx = np.repeat(inout, ngals_tot)
            ind_ok = np.where(lambda_evol_all[:,inout] > 0.01)[0]

            if len(ind_ok) > 0:
                yy = lambda_evol_all[ind_ok,inout]#[:,np.newaxis]
                #print(inout, len(yy))
                #xx,yy,zz = density_map_slow(xx[ind_ok], all_data[ind_ok])
                zz = gaussian_kde(yy)(yy)
                #kde = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(X)
                #log_dens = kde.score_samples(X_plot)
                zz /= max(zz)

                idx = zz.argsort()
                #yy = yy[idx]
                #zz = zz[idx]

                all_xxs.extend(xx[idx])
                all_yys.extend(yy[idx])
                all_zzs.extend(zz[idx])

                #patches.append(Circle((xx[ind_ok], yy), 60,
                #                color=cmap.zz,
                #                facecolor='none'))
                                #edgecolor=(0, 0.8, 0.8),
                                #linewidth=3, alpha=0.5))

                #ax.scatter(xx[ind_ok], yy, c=zz, s=60,
                #          edgecolor='',
                #          cmap=cmap,
                #          rasterized=True)

            print(nout, end='\r')
            ax.scatter(all_xxs, all_yys, c=all_zzs, s=80,
                  edgecolor='',
                  cmap=cmap,
                  rasterized=True)

        #patches.
        #p = PatchCollection(patches)
        #p.set_array(np.array(colors))
        #ax.add_collection()

        lambda_range=[0.01, 0.8]
        yticks_ok=[0.0, 0.2, 0.4, 0.6, 0.8]
        ax.set_ylim([-0.01, 0.9])
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([str(yy) for yy in yticks_ok])

    ax.set_xlim([nnouts+1,-1])

    #plt.show()
    if fname is not None:
        # Common, but if ax is to be returned,
        # I will want to have control over labels. (for subplots)
        ax.set_ylabel(r"$\lambda_{R_{eff}}$", fontsize=20)
        ax.set_xlabel("Redshift")

        plt.tight_layout()
        plt.savefig(fname, dpi=200, rasterized=True)
    else:
        print("\n FNAME", fname)

    #return im
