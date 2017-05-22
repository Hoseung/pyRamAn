
# coding: utf-8

# In[1]:

import load
import numpy as np
import matplotlib.pyplot as plt
from galaxy import galaxy
import tree


wdir = './29172/'


galcat = tree.halomodule.Halo(base=wdir, is_gal=True, nout=187)
galcat.data["id"][galcat.data["m"] > 5e10]


nout=187
gind = 10
gal = galaxy.Galaxy(halo=galcat.data[gind], info=galcat.info)
gm = load.rd_GM.rd_gal(nout, galcat.data["id"][gind], wdir=wdir)


#print("")
print(len(gm.star))
gal.mk_gal(gm.star,None,None,unit_conversion="GM")
#print(len(gal.star))


lambdas = gal.cal_lambda_r_eps(save_result=False)

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

plt.imshow(gal.vmap, origin="lower")
plt.scatter(xinds[i_ok], yinds[i_ok], color="g")
plt.show()

v_good = gal.vmap.flatten()[i_ok]
sig_good = gal.sigmap.flatten()[i_ok]
d_good = np.sqrt(np.square(xinds[i_ok]-fit["xcen"]) + np.square(yinds[i_ok]-fit["ycen"]))
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


star = gal.star
i_inner = np.where((np.square(star['x']) +
                    np.square(star['y']) +
                    np.square(star['z'])) < gal.meta.reff**2)[0]# in kpc unit

sig = weighted_std(star["vz"], star['m'])

binsize=5
v_smooth=np.array([np.mean(v_good[i*binsize:(i+1)*binsize]) for i in range(np.int(len(v_good)/binsize))])
d_smooth=np.array([np.mean(d_good[i*binsize:(i+1)*binsize]) for i in range(np.int(len(v_good)/binsize))])

from scipy.signal import argrelmax
from scipy.interpolate import UnivariateSpline


d_sort = np.argsort(d_good)
v_good = np.abs(v_good[d_sort])
sig_good= sig_good[d_sort]
d_good = d_good[d_sort]

plt.scatter(d_good, v_good, color="r")
plt.scatter(d_good, sig_good, color="b")

imax = argrelmax(v_smooth)[0][0]
vmax = v_smooth[imax]
plt.plot(d_smooth, v_smooth, color="r")
plt.scatter(d_smooth[imax], v_smooth[imax], s=200, marker="^", color="r")
plt.scatter(0, sig, s=200, marker="v", color='b')

print(vmax, sig, vmax/sig)
plt.show()


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


star = gal.star
i_inner = np.where((np.square(star['x']) +
                    np.square(star['y']) +
                    np.square(star['z'])) < gal.meta.reff**2)[0]# in kpc unit

sig = weighted_std(star["vz"], star['m'])

binsize=5
v_smooth=np.array([np.mean(v_good[i*binsize:(i+1)*binsize]) for i in range(np.int(len(v_good)/binsize))])
d_smooth=np.array([np.mean(d_good[i*binsize:(i+1)*binsize]) for i in range(np.int(len(v_good)/binsize))])

from scipy.signal import argrelmax
from scipy.interpolate import UnivariateSpline

d_sort = np.argsort(d_good)
v_good = np.abs(v_good[d_sort])
sig_good= sig_good[d_sort]
d_good = d_good[d_sort]

plt.scatter(d_good, v_good, color="r")
plt.scatter(d_good, sig_good, color="b")
#v_smooth = smooth(v_good, beta=10, window_len=50)
#v_smooth = UnivariateSpline(d_good, v_good, k=5)
#sig_smooth = smooth(sig_good, beta=10, window_len=30)

imax = argrelmax(v_smooth)[0][0]
vmax = v_smooth[imax]
#plt.plot(d_good, v_smooth(d_good), color="r")
plt.plot(d_smooth, v_smooth, color="r")
plt.scatter(d_smooth[imax], v_smooth[imax], s=200, marker="^", color="r")
plt.scatter(0, sig, s=200, marker="v", color='b')

#plt.plot(d_good, sig_smooth, color="b")
print(vmax, sig, vmax/sig)
plt.savefig()

          
