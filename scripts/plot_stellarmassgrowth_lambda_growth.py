# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 02:57:47 2015

1.
This script plots co-evolution history of galaxy stellar mass growth and 
lambda_r evolution. Considering abrupt stellar mass growth as a consequence of 
galaxy merger, the role of merger history on lambda_r might be understood.

2.
This script loads lambda_mp.py output catalogs (.pickle).

3. 
For the moment, only galaxies with full halo merger history are stored in catalogs. 

2015/08/18

Multiple yaxis 

@author: hoseung
"""

import pickle
import numpy as np

wdir = '/home/hoseung/Work/data/05427/'

nout_fi = 187
nout_ini = 37
nouts = np.arange(nout_ini, nout_fi + 1)
nnouts = len(nouts)



def load_pickle(fname):
    with open(fname, 'rb') as f:
        return pickle.load(f)

# Final galaxies
# Why set..?
#final_gals = set()
for inout, nout in enumerate(reversed(np.arange(nout_ini, nout_fi + 1))):
    cat = load_pickle(wdir + 'catalog/' + 'catalog' + str(nout) + '.pickle')
    #final_gals.update(cat['final'])
    final_gals = cat['final_gal']

ngals = len(final_gals)
mstar = np.zeros((ngals, nnouts))
l_r = np.zeros((ngals, nnouts))
reff = np.zeros((ngals, nnouts))
fg = np.zeros((ngals, nnouts), dtype=int)
print(ngals)
final_gals = np.sort(list(final_gals))

#%%#######################################################################
# Read catalogs and extract mstar and lambda-r of galaxies at all nouts.
for inout, nout in enumerate(reversed(np.arange(nout_ini, nout_fi + 1))):
    cat = load_pickle(wdir + 'catalog/' + 'catalog' + str(nout) + '.pickle')
    for igal, idgal in enumerate(cat['final_gal']):
        #ind = which galaxy..? -> need tree. (save the ID of final halo too)           
        ind = np.where(idgal == final_gals)[0]
#        print(ind, final_gals[ind], inout, nout, idgal, cat['id'][igal])
        if len(ind) > 0 :
            fg[ind,inout] = final_gals[ind]
            mstar[ind,inout] = cat['mstar'][igal]
            l_r[ind,inout] = cat['lambda_r'][igal]
            reff[ind,inout] = cat['rgal'][igal]


#%%
zreds=[]
aexps=[]
import load
for nout in nouts:
    info = load.info.Info(nout=nout, base=wdir, load=True)
    aexps.append(info.aexp)
    zreds.append(info.zred)
aexps = np.array(aexps)
zreds = np.array(zreds)

#%%
def aexp2zred(aexp):
    return [1.0/a - 1.0 for a in aexp]

def zred2aexp(zred):
    return [1.0/(1.0 + z) for z in zred]

def lbt2aexp(lts):
    import astropy.units as u
    from astropy.cosmology import WMAP7, z_at_value
    zreds = [z_at_value(WMAP7.lookback_time, ll * u.Gyr) for ll in lts]
    return [1.0/(1+z) for z in zreds]

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
mm = mstar/1e10

plt.close()
plt.ioff()

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

for i, idgal in enumerate(cat['final_gal']):
#for i, idgal in enumerate([1618]):
    if mm[i][0] < 0.2 :
        print(idgal, mm[i][0], mm[i][-1])
        continue
    print(idgal, "!!!!")

    plt.rcParams["figure.figsize"] = [12,10]
    fig, axes = plt.subplots(3)
#    plt.figure(num=1, figsize=[10,20])
    fig.suptitle("ID: " + str(idgal).zfill(5), fontsize=18)#, y=1.01)
    lns1 = axes[0].plot(nouts[::-1], l_r[i], label=r"$\lambda_{R}$")
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
    lns2 = axes[1].plot(nouts[::-1], mm[i], 'r-', label="stellar mass")
    axes[1].set_ylim([0, 1.3*max(mm[i])])
    axes[1].set_xlim([37,187])
    axes[1].set_ylabel(r"Stellar mass $[10^{10}M_{\odot}]$")
    axes[1].get_yaxis().get_offset_text().set_y(1)
    
#    ax3 = ax1.twinx() # Reff
#    ax3.spines["right"].set_position(("axes", 1.2))
#    make_patch_spines_invisible(ax3)
    # Second, show the right spine.
#    ax3.spines["right"].set_visible(True)
    axes[2].set_ylabel("Reff [kpc]")
    axes[2].set_xlim([37,187])
    lns3 = axes[2].plot(nouts[::-1], reff[i], 'g-', label='Reff')

    # hide x axes so that subplots stick together.
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    
    ax4 = axes[0].twiny()
    ax4.set_xlabel("Lookback time", labelpad=10)
    ax4.set_xticks(lbt_pos)
    ax4.set_xticklabels(lbt_target_str)
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    axes[0].legend(lns, labs, loc=0)
    # logend location codes:
    # 0 ~ 10 
    # best, ur, ul, lr, ll, r, cl, cr, lower c, upper c, center
    #
       
#    plt.show()
    plt.savefig(wdir + 'catalog/' + str(idgal).zfill(5) + '.png')
    plt.close()