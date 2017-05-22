
# coding: utf-8

# In[1]:

import load
import numpy as np
import matplotlib.pyplot as plt
from galaxymodule import galaxy
import tree
from scipy.signal import argrelmax
from scipy.interpolate import UnivariateSpline

s = load.sim.Sim(nout=782)
s.add_part()


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


wdir = './'

nout=315
galcat = tree.halomodule.Halo(base=wdir, is_gal=True, nout=nout)

vmaxs=[]
sigs=[]
ratios=[]
ids=[]
lambda_rs=[]
gals=[]
n_pseudo = 1
nreff=1.0
#masses=[]

large_galaxies = galcat.data["id"][galcat.data["m"] > 1e9]
#large_galaxies = np.intersect1d(large_galaxies, np.arange(5,124736,10))

print(len(large_galaxies))

npix_per_reff=5

for gid in large_galaxies:
    gind = np.where(galcat.data["id"] == gid)[0]

    gal = galaxy.Galaxy(halo=galcat.data[gind], info=galcat.info)
    gm = load.rd_GM.rd_gal(nout, galcat.data["id"][gind][0], wdir=wdir)

    #print("")
    if len(gm.star) < 2e3:
        print("skip")
        continue
    good_gal = gal.mk_gal(gm.star,None,None,unit_conversion="GM",
                           den_lim=1e6, den_lim2=5e6)
    #                      den_lim=1e8, den_lim2=5e8)
    if not good_gal or gal.meta.mstar < 5e9:
        continue
    #print(len(gal.star))

    lambdas = gal.cal_lambda_r_eps(save_result=False, n_pseudo=n_pseudo,
                                   npix_per_reff=npix_per_reff)
    fit = gal.meta.mge_result_list[0] # Use ellipse fitting at 1Reff, because most of the time Vm is inside 1Reff.

    # get points near the major axis line
    img_size = 40
    xinds = np.tile(np.arange(img_size) + 0.5, img_size)
    yinds = np.repeat(np.arange(img_size) + 0.5, img_size)

    # polynomial constants
    f_a = np.tan(fit["pa_rad"])
    f_b = -1
    f_c = fit["ycen"] - f_a*fit["xcen"]

    distance_to_line = np.abs(f_a*xinds + f_b*yinds + f_c)/np.sqrt(f_a**2 + f_b**2)
    i_ok = distance_to_line < 1

    fig, ax = plt.subplots(2)

    ax[0].imshow(gal.vmap, origin="lower")
    ax[0].scatter(xinds[i_ok], yinds[i_ok], color="g")
    #plt.savefig(str(gid) + "_vmap.png")
    #plt.show()

    v_good = gal.vmap.flatten()[i_ok]
    sig_good = gal.sigmap.flatten()[i_ok]
    d_good = np.sqrt(np.square(xinds[i_ok]-fit["xcen"]) + np.square(yinds[i_ok]-fit["ycen"]))


    d_sort = np.argsort(d_good)
    v_good = np.abs(v_good[d_sort])
    sig_good= sig_good[d_sort]
    d_good = d_good[d_sort]

    plt.scatter(d_good, v_good, color="r")
    plt.scatter(d_good, sig_good, color="b")

    npoints = len(v_good)
    binsize = np.int(npoints/10)

    # smoothed arrays
    v_smooth=np.array([np.mean(v_good[i*binsize:(i+1)*binsize]) for i in range(10)])
    d_smooth=np.array([np.mean(d_good[i*binsize:(i+1)*binsize]) for i in range(10)])

    #i_vmax = argrelmax(v_smooth)[0][0]
    #r_vmax = d_smooth[i_vmax]
    #r_at_max = min([r_vmax,2*gal.meta.reff])
    #print(d_smooth)
    #print(nreff, gal.meta.reff)
    #print
    imax = np.argmax(v_smooth[d_smooth<nreff*npix_per_reff])
    # img_size = rgal/reff
    # maximum velocity at ...?
    #imax = np.argmax(d_smooth>r_at_max)
#    print(imax)
#    print("_______")

    # sigma
    star = gal.star
    star = gal.star[np.where((np.square(star['x']) +
                        np.square(star['y']) +
                        np.square(star['z'])) < gal.meta.reff**2)[0]]# in kpc unit
    #star = star[i_inner]

    sig = weighted_std(star["vz"], star['m'])

    vmax = v_smooth[imax]
    plt.plot(d_smooth, v_smooth, color="r")
    plt.scatter(d_smooth[imax], v_smooth[imax], s=200, marker="^", color="r")
    plt.scatter(0, sig, s=200, marker="v", color='b')

    print(vmax, sig, vmax/sig)
    ax[1].set_aspect('auto')
    plt.savefig(str(gid) + "_vel_curve" + str(n_pseudo) + ".png")

    gal.meta.vmax = vmax
    gal.meta.sig = sig
    gal.meta.v_sig = vmax/sig

    vmaxs.append(vmax)
    sigs.append(sig)
    ratios.append(vmax/sig)
    ids.append(gid)
    lambda_rs.append(gal.meta.lambda_r[0])
    gals.append(gal)
    plt.close()


import pickle

ratios = np.array(ratios)
lambda_rs = np.array(lambda_rs)

pickle.dump((ids,vmaxs,sigs,ratios, lambda_rs), open("results_" + str(n_pseudo) + ".pickle","wb"))
pickle.dump(gals, open("gals_results_" + str(n_pseudo) + ".pickle","wb"))

i_ok = lambda_rs < 0.9

fig, ax = plt.subplots()

ax.scatter(ratios[i_ok], lambda_rs[i_ok])
ax.set_ylabel(r"$\lambda$")
ax.set_xlabel(r"$v/\sigma$")
ax.set_xlim([0,1.3])
ax.set_ylim([0,0.9])
plt.savefig("Lambda_vs_vsig_" + str(n_pseudo) + ".png")


fig, ax = plt.subplots(1,2)
fig.set_size_inches(8,4)

ax[0].scatter(np.log10([gal.meta.mstar for gal in gals]), ratios[i_ok])
ax[1].scatter(np.log10([gal.meta.mstar for gal in gals]), lambda_rs[i_ok])
ax[0].set_ylabel(r"$v/\sigma$")
ax[1].set_ylabel(r"$\lambda$")

ax[0].set_xlabel(r"$log(M_{*})$")
ax[1].set_xlabel(r"$log(M_{*})$")
#ax.set_xlim([0,1.3])
#ax.set_ylim([0,0.9])
plt.tight_layout()
plt.savefig("Lambda_vsig_mass_" + str(n_pseudo) + ".png")

print("Number of galaxies", len(gals))
