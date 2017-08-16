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
    zz_target = [0, 0.5, 1.0, 2.0, 3.0]

    x_tick_pos = match.match_list_ind(nouts, np.array([782, 535, 343, 197, 125]))# - nouts[-1]

    ax.set_xlabel("Redshift")
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
    zz_target = [0, 0.5, 1.0, 2.0, 3.0]

    x_tick_pos = match.match_list_ind(nouts, np.array([782, 535, 343, 197, 125]))

    nnouts = len(nouts)

    # For a given list of nouts,
    # calculate a nice-looking set of zreds AND lookback times
    x_tick_pos = np.searchsorted(zreds[::-1], zz_target)[::-1]# + nout_ini# + nout_min
    # searchsorted requires arrays be in ascending order.

    #  at z = 3, lbt = 11.5243, age of universe = 2.18844
    u_age_targets=[3,5,8,10,13]
    z_targets_u_age = age2zred(u_age_targets)
    u_age_pos = np.searchsorted(zreds, z_targets_u_age)

    ax2.set_xticks(u_age_pos)
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
                     nnza,
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


    # compile data
    nnouts = len(nouts)
    ngals_tot = len(serial_results)

    lambda_evol_all = np.zeros([ngals_tot, nnouts])
    nstep_max = nnza.nnza["nstep"].max()
    # Starting nouts of galaxies are different.
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

    fig, ax = plt.subplots(1)
    ax2 = ax.twiny()

    modify_ticks1(ax, nnza, nouts)
    modify_ticks2(ax2, nnza, nouts)

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
                       gridsize=50,
                       bins=None,
                       cmap=cmap)
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
        ax.set_ylabel(r"$\lambda_{R_{eff}}$", fontsize=20)

    #plt.show()
    if fname is not None:
        plt.savefig(fname)

    return im
