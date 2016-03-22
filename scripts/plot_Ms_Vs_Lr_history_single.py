# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 09:08:08 2015

sinlge galaxy vversion of plot_stellarmassgrowth_lambda_growth.py

@author: hoseung
"""
import pickle
import numpy as np

def load_pickle(fname):
    with open(fname, 'rb') as f:
        return pickle.load(f)


wdir = '/home/hoseung/Work/data/05427/'

nout_fi = 187
nout_ini = 37
nouts = np.arange(nout_ini, nout_fi + 1)
nnouts = len(nouts)
galid = 1392
# Final galaxies
# Why set..?
#final_gals = set()
cat = load_pickle(wdir + 'catalog/individual/catalog187_' + str(galid) + '.pickle')
# Starts from 187 and goes backwards.

def aexp2zred(aexp):
    return [1.0/a - 1.0 for a in aexp]

def zred2aexp(zred):
    return [1.0/(1.0 + z) for z in zred]

def lbt2aexp(lts):
    import astropy.units as u
    from astropy.cosmology import WMAP7, z_at_value
    zreds = [z_at_value(WMAP7.lookback_time, ll * u.Gyr) for ll in lts]
    return [1.0/(1+z) for z in zreds]

zreds=[]
aexps=[]
import load
for nout in nouts:
    info = load.info.Info(nout=nout, base=wdir, load=True)
    aexps.append(info.aexp)
    zreds.append(info.zred)
aexps = np.array(aexps)
zreds = np.array(zreds)


# For a given list of nouts, 
# calculate a nice-looking set of zreds.
# AND lookback times
z_targets=[0, 0.2, 0.5, 1, 2, 3]
z_target_str=["{:.2f}".format(z) for z in z_targets]
a_targets_z = zred2aexp(z_targets)
z_pos =  [nout_ini + (1 - (max(aexps) - a)/aexps.ptp()) * nnouts for a in a_targets_z]

lbt_targets=[0.00001,1,3,5,8,12]
lbt_target_str=["{:.0f}".format(l) for l in lbt_targets]
a_targets_lbt = lbt2aexp(lbt_targets)
lbt_pos = [nout_ini + (1 - (max(aexps) - a)/aexps.ptp()) * nnouts for a in a_targets_lbt]
#from astropy.cosmology import WMAP7 as cosmo
#lookback_t=[cosmo.lookback_time(i).value for i in zreds]

#%%
import matplotlib.pyplot as plt
# plot each galaxy.
# stellar mass growth and lambda_r as a function of time.

# The exponent (also called ass "offset") in the figure (1e11)
# overlaps with lookback time tick labels.
# And moving the offset around is not easy. 
# So, manually divide the values.
#mm = mstar/1e10

#plt.close()
#plt.ioff()

def make_patch_spines_invisible(ax):
    """
        Useful for plotting multiple variables (more than two twinx())
    """
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    #for sp in ax.spines.itervalues():  
    # Changed in Python3
    for sp in ax.spines.values():
        sp.set_visible(False)


l_r =[]
for i in cat['lambda_arr']:
    l_r.append(i[4])

mm = cat['mstar']/1e10
#l_r = cat['lambda_r']
reff = cat['rgal']

plt.rcParams["figure.figsize"] = [12,10]
fig, axes = plt.subplots(2)

#    plt.figure(num=1, figsize=[10,20])
fig.suptitle("ID: " + str(galid).zfill(5), fontsize=18)#, y=1.01)
nn = nouts#[::-1]
lns1 = axes[0].plot(nn[1:], l_r, label=r"$\lambda_{R}$")
axes[0].set_xticks(z_pos)
axes[0].set_xticklabels(z_target_str)
plt.subplots_adjust(left = 0.1, right = 0.9, \
                    wspace = 0.1, hspace = 0.0, \
                    bottom = 0.1, top = 0.85)

axes[0].set_xlim([37,187])
axes[0].set_ylim([0,1.0])
axes[0].set_ylabel(r"$\lambda_{R}$")
axes[0].set_xlabel("redshift")

#    ax2 = axes[0].twinx()
lns2 = axes[1].plot(nn[1:], mm, 'r-', label="stellar mass")
axes[1].set_ylim([0, 1.3*max(mm)])
axes[1].set_xlim([37,187])
axes[1].set_ylabel(r"Stellar mass $[10^{10}M_{\odot}]$")
axes[1].get_yaxis().get_offset_text().set_y(1)

#    ax3 = ax1.twinx() # Reff
#    ax3.spines["right"].set_position(("axes", 1.2))
#    make_patch_spines_invisible(ax3)
# Second, show the right spine.
#    ax3.spines["right"].set_visible(True)
"""
axes[2].set_ylabel("Reff [kpc]")
axes[2].set_xlim([37,187])


lns3 = axes[2].plot(nn[1:], reff, 'g-', label='Reff')

# hide x axes so that subplots stick together.
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
"""
ax4 = axes[0].twiny()
ax4.set_xlabel("Lookback time", labelpad=10)
ax4.set_xticks(lbt_pos)
ax4.set_xticklabels(lbt_target_str)
lns = lns1+lns2#+lns3
labs = [l.get_label() for l in lns]
axes[0].legend(lns, labs, loc=0)
# logend location codes:
# 0 ~ 10 
# best, ur, ul, lr, ll, r, cl, cr, lower c, upper c, center
#
   
plt.show()
    #plt.savefig(wdir + 'catalog/' + str(idgal).zfill(5) + '.png')
    #plt.close()
