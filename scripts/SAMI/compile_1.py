# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 13:30:38 2015

@author: hoseung
"""

import pickle
import numpy as np

wdir = '/home/hoseung/Work/data/05427/'

nout_fi = 187
nout_ini = 120
nouts = np.arange(nout_ini, nout_fi + 1)
nnouts = len(nouts)

# Final galaxies
final_gals = set()
for inout, nout in enumerate(reversed(np.arange(nout_ini, nout_fi + 1))):
    with open(wdir + 'catalog/' + 'catalog' + str(nout) + '.pickle', 'rb') as f:
        cat = pickle.load(f)
        final_gals.update(cat['final_gal'])

ngals = len(final_gals)
mstar = np.zeros((ngals, nnouts))
l_r = np.zeros((ngals, nnouts))
print(ngals)

final_gals = np.sort(list(final_gals))
#%%
for inout, nout in enumerate(reversed(np.arange(nout_ini, nout_fi + 1))):
    with open(wdir + 'catalog/' + 'catalog' + str(nout) + '.pickle', 'rb') as f:
        cat = pickle.load(f)
        for igal, idgal in enumerate(cat['final_gal']):
            #ind = which galaxy..? -> need tree. (save the ID of final halo too)           
            ind = np.where(idgal == final_gals)[0]
            if len(ind) > 0 :
                mstar[ind,inout] = cat['mstar'][igal]
                l_r[ind,inout] = cat['lambda_r'][igal]

#%%    

# plot each galaxy.
# stellar mass growth and lambda_r as a function of time.
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

# For a given list of nouts, 
# calculate a nice-looking set of zreds.
# AND lookback times
z_targets=[0, 0.2, 0.5, 1, 2, 3]
z_target_str=["{:.2f}".format(z) for z in z_targets]

a_targets_z = zred2aexp(z_targets)
z_pos =  [nout_ini + (1 - (max(aexps) - a)/aexps.ptp()) * nnouts for a in a_targets_z]
#%%
def lbt2aexp(lts):
    import astropy.units as u
    from astropy.cosmology import WMAP7, z_at_value
    zreds = [z_at_value(WMAP7.lookback_time, ll * u.Gyr) for ll in lts]
    return [1.0/(1+z) for z in zreds]

lbt_targets=[0.00001,1,3,5,8,12]
lbt_target_str=["{:.0f}".format(l) for l in lbt_targets]

a_targets_lbt = lbt2aexp(lbt_targets)
lbt_pos = [nout_ini + (1 - (max(aexps) - a)/aexps.ptp()) * nnouts for a in a_targets_lbt]
#from astropy.cosmology import WMAP7 as cosmo
#lookback_t=[cosmo.lookback_time(i).value for i in zreds]
#Ipr
#

#%%
#z_pos=[113, 88, 63, 38, 16, 0]
#z_pos = [0.1, 0.2, 0.4, 0.5, 0.6, 0.9]

# The exponent (also called ass "offset") in the figure (1e11)
# overlaps with lookback time tick labels.
# And moving the offset around is not easy. 
# So, manually divide the values.
mm = mstar/1e10

import matplotlib.pyplot as plt

plt.close()
plt.ioff()


fig, ax1 = plt.subplots()
fig.subplots_adjust(right=0.75)

ax1.set_xticks(z_pos)
ax1.set_xticklabels(z_target_str)

ax1.set_xlim([37,187])
ax1.set_ylim([0,1.0])
ax1.set_ylabel(r"$\lambda_{R}$")
ax1.set_xlabel("redshift")

ax2 = ax1.twinx()    
ax2.set_ylim([0, 10.0])
ax2.set_ylabel(r"Stellar mass $[10^{10}M_{\odot}]$")
ax2.get_yaxis().get_offset_text().set_y(1)

ax3 = ax1.twiny()
ax3.set_xlabel("Lookback time")
ax3.set_xticks(lbt_pos)
ax3.set_xticklabels(lbt_target_str)
#lns = lns1+lns2
#labs = [l.get_label() for l in lns]
#ax1.legend(lns, labs, loc=0)
ax1.plot(nouts, l_r[0], label=r"$\lambda_{R}$")
m_final = mm[0][-1]
ax2.plot(nouts[::-1], mm[0]/m_final, 'r-', label="stellar mass")

for i, idgal in enumerate(cat['final_gal'][10:11]):
    m_final = mm[i][-1]
    print(m_final)
    if m_final < 1 :
        continue
        #lns1 = ax1.plot(nouts[::-1], l_r[i], label=r"$\lambda_{R}$")    
    # ratio. 
    mthis=mm[i]/m_final
    print(mthis)
    #ax2.plot(nouts, mm[i]/m_final, 'r-', label="stellar mass")
    ax2.scatter(nouts[::-1], mm[i]/m_final, label="stellar mass")
    
plt.show()
#    plt.savefig(wdir + 'catalog/' + str(idgal).zfill(5) + '.png')
#    plt.close()