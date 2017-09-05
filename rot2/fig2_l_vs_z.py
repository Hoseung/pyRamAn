## Fig2
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP7, z_at_value
from utils import match

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

    aexps = nnza.a2b(nouts, "nout", "aexp")
    zreds = nnza.a2b(nouts, "nout", "zred")
    nnouts = len(nouts)

    # For a given list of nouts,
    # calculate a nice-looking set of zreds AND lookback times
    zz_target = [0.0, 0.5, 1.0, 2.0]
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
    aexps = nnza.a2b(nouts, "nout", "aexp")
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

def density_map(x, y, sort=True):
    from scipy.stats import gaussian_kde
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    z /= max(z)

    idx = z.argsort()
    xx, yy = x[idx], y[idx]
    z = z[idx]

    #im = ax.scatter(xx, yy, c=z, s=50, edgecolor='')
    return xx,yy,z


def plot_lambda_evol(serial_results, nouts,
                     nnza, ax=None,
                     density = "hexbin",
                     fname=None,
                     cmap ="jet"):
    """
    im= plot_lambda_evol(serial_results, nouts,
                 nnza_cell,
                 density = "kernel",
                 fname="./RUN2/figs/Lambda_evol.png",
                 cmap ="jet")
    """

    if hasattr(serial_results[0], "finearr"):
        fine=True
        #print("fine")
    else:
        fine=False

    # compile data
    nnouts = len(nouts)
    ngals_tot = len(serial_results)

    lambda_evol_all = np.zeros([ngals_tot, nnouts])
    nstep_max = nnza.nnza["nstep"].max()
    nstep_min = nstep_max - nnouts
    # Starting nouts of galaxies are different.
    if fine:
        for igal, gal in enumerate(serial_results):
            if len(gal.finearr) == 0:
                continue
            if True:
                nstep = gal.finearr["nstep"]
                #print(igal, nstep_max-nstep, nstep_min)
                ind = np.where(gal.finearr["nstep"] > nstep_min)[0]
                lambda_evol_all[igal,ind] = gal.finearr["lambda_r"][ind]
            else:
                pass
    else:
        for igal, gal in enumerate(serial_results):
            if len(gal.data) == 0:
                continue
            for gg in gal.data[0][0]:
                try:
                    nstep = gg.nstep
                    lambda_evol_all[igal][nstep_max - nstep] = gg.lambda_r[0]
                    #print("Good")
                except:
                    pass

    # change ticks
    aexps = nnza.a2b(nouts, "nout", "aexp")
    zreds = nnza.a2b(nouts, "nout", "zred")

    if ax is None:
        fig, ax = plt.subplots()
        fig.set_size_inches(6,4)
        fname=None
    #ax2 = ax.twiny()

    modify_ticks1(ax, nnza, nouts)
    #modify_ticks2(ax2, nnza, nouts)

    lambda_range = nouts.ptp()

    if density == "heat":
        den_map = np.zeros((nbins, nnouts))
        for i in range(nnouts):
            den_map[:,i], ypoints = np.histogram(lambda_evol_all[:,i], bins=nbins, range=lambda_range)
            den_map[:,i] /= den_map[:,i].max()
        im = ax.imshow(den_map, origin="lower",
                       cmap=cmap,
                       interpolation="gaussian")
    elif density == "hexbin":
        xx = np.tile(np.arange(nnouts), ngals_tot)
        all_data = lambda_evol_all.ravel()
        ind_ok = all_data > 0.01
        #xx,yy,z = density_map(xx[ind_ok], all_data[ind_ok])
        im = ax.hexbin(xx[ind_ok], 10 * all_data[ind_ok],
                       gridsize=40,
                       bins=None,
                       cmap=cmap)

        y_tick_pos=np.arange(10)
        ax.set_yticks(y_tick_pos)#[::-1]
        ax.set_yticklabels(labels = ["{:0.1f}".format(0.1*y) for y in y_tick_pos])


    elif density == "kernel":
        xx = np.tile(np.arange(nnouts), ngals_tot)
        all_data = lambda_evol_all.ravel()
        ind_ok = all_data > 0.01
        xx,yy,z = density_map(xx[ind_ok], all_data[ind_ok])
        im = ax.scatter(xx, yy, c=z, s=50, edgecolor='', cmap=cmap)
        lambda_range=[0.01, 0.8]
        yticks_ok=[0.0, 0.2, 0.4, 0.6, 0.8]
        ax.set_ylim([-0.05, 0.9])
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
        plt.savefig(fname)

    return im
