# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 14:20:02 2015
Plot Delta Lambda_r Vs Delta Stellar mass


@author: hoseung
"""

import matplotlib.pyplot as plt
import pickle
import numpy as np

## fucntions 
def load_pickle(fname):
    with open(fname, 'rb') as f:
        return pickle.load(f)

## time
def aexp2zred(aexp):
    return [1.0/a - 1.0 for a in aexp]

def zred2aexp(zred):
    return [1.0/(1.0 + z) for z in zred]

def lbt2aexp(lts):
    import astropy.units as u
    from astropy.cosmology import WMAP7, z_at_value
    zreds = [z_at_value(WMAP7.lookback_time, ll * u.Gyr) for ll in lts]
    return [1.0/(1+z) for z in zreds]

def density_map(x, y, ax, sort=True):
    from scipy.stats import gaussian_kde
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy) 
    z /= max(z)

    idx = z.argsort()    
    xx, yy = x[idx], y[idx]
    z = z[idx]
    
    im = ax.scatter(xx, yy, c=z, s=50, edgecolor='')
    return im

def smooth(x,beta):
 """ kaiser window smoothing """
 window_len=11
 # extending the data at beginning and at the end
 # to apply the window at the borders
 s = np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
 w = np.kaiser(window_len,beta)
 y = np.convolve(w/w.sum(),s,mode='valid')
 return y[5:len(y)-5]

def outlier_clip(data, tol, max_size_bad = 10):
    for dx in range(2*int(max_size_bad/2)-1, 8, -2):
        # first test must be with dx = 1
        #print(dx)
        ll = data[:-2*dx]
        cc = data[dx:-dx]
        rr = data[2*dx:]
       
        bad1 = np.logical_and( (cc-ll)/cc > tol, (cc-rr)/cc > tol )
        bad2 = np.logical_and( (cc-ll)/cc < -tol, (cc-rr)/cc < -tol )
        bad = np.logical_or(bad1 , bad2)
        
        cc[bad] = 0.5 * (ll[bad] + rr[bad])    

    return data


def fill_inf(data):
    # fill infinite element with the next element
    ind = np.isinf(data)
    data[ind] = data[1:][ind]
    return data


#%%
# multiple clusters, but with fixed nout range.

nout_fi = 187
nout_ini = 150
nouts = np.arange(nout_ini, nout_fi + 1)
nnouts = len(nouts)

clusters = [5427, 36413, 74010, 1605][0:1]
exclude_gals = [[1496, 85, 1636, 1340],[],[],[]]

all_l_r =np.zeros(0)
all_zreds = np.zeros(0)
cluster_data=[]
for i, cluster in enumerate(clusters):
    wdir = '/home/hoseung/Work/data/' + str(cluster).zfill(5) 

    if i == 0:
        ##  calculate time
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
        # calculate a nice-looking set of zreds AND lookback times
        z_targets=[0, 0.2, 0.5, 1, 2, 3]
        z_target_str=["{:.2f}".format(z) for z in z_targets]
        a_targets_z = zred2aexp(z_targets)
        z_pos =  [nout_ini + (1 - (max(aexps) - a)/aexps.ptp()) * nnouts for a in a_targets_z]
        
        lbt_targets=[0.00001,1,3,5,8,12]
        lbt_target_str=["{:.0f}".format(l) for l in lbt_targets]
        a_targets_lbt = lbt2aexp(lbt_targets)
        lbt_pos = [nout_ini + (1 - (max(aexps) - a)/aexps.ptp()) * nnouts for a in a_targets_lbt]
        
        #lookback_t=[cosmo.lookback_time(i).value for i in zreds]


    # Load catalog
    #for inout, nout in enumerate(reversed(np.arange(nout_ini, nout_fi + 1))):
    
    cat = load_pickle(wdir + '/catalog/' + 'catalog' + str(nout_fi) + '.pickle')
    final_gals = list(cat['final_gal'])
    
    # exclude disky galaxies
    for bad_gal in exclude_gals[i]:
        final_gals.remove(bad_gal)
    
    ngals = len(final_gals)
    mstar = np.zeros((ngals, nnouts))
    l_r = np.zeros((ngals, nnouts))
    reff = np.zeros((ngals, nnouts))
    time = np.zeros((ngals, nnouts))
    fg = np.zeros((ngals, nnouts), dtype=int)
    
    final_gals = np.sort(list(final_gals))
    
    # Read catalogs and extract mstar and lambda-r of galaxies at all nouts.
    #    for inout, nout in enumerate(reversed(np.arange(nout_ini, nout_fi + 1))):
    for inout, nout in enumerate(np.arange(nout_ini, nout_fi + 1)):
#        print(nout)
        cat = load_pickle(wdir + '/catalog/' + 'catalog' + str(nout) + '.pickle')
        
        for igal, idgal in enumerate(cat['final_gal']):
            #ind = which galaxy..? -> need tree. (save the ID of final halo too)           
            ind = np.where(idgal == final_gals)[0]
            
            if len(ind) > 0 :
                fg[ind,inout] = final_gals[ind]
                mstar[ind,inout] = cat['mstar'][igal]
                l_r[ind,inout] = cat['lambda_r'][igal]
                reff[ind,inout] = cat['rgal'][igal]
                time[ind,inout] = zreds[inout]

    cluster_data.append({})
    cluster_data[i].update({"ngals":ngals, "mstar": mstar, "l_r":l_r,
                        "reff":reff,"time":time,"fg":fg, "final_gals":final_gals})

    #zzz = np.repeat(zreds, ngals)
    zzz = time.ravel()
    all_zreds = np.concatenate((all_zreds, zzz))
    #print(all_aexps, zzz)
    all_l_r = np.concatenate((all_l_r, l_r.ravel()))#, axis=1)

from astropy.cosmology import WMAP7 as cosmo
lookback_t=np.asarray([cosmo.lookback_time(i).value for i in all_zreds])


#%%
# smooth
# [:,0] = 187, [:,1] = 186
dt = 3
dx_m = np.zeros((ngals,nnouts - dt))
dx_l = np.zeros((ngals,nnouts - dt))
arrnouts = np.zeros((ngals,nnouts -dt))
for i in range(nnouts-dt):
    dx_m[:,i] = (mstar[:,i+dt] - mstar[:,i])/mstar[:,i]
    dx_l[:,i] = l_r[:,i+dt] - l_r[:,i]
    arrnouts[:][i] = i

# plot

fig, ax = plt.subplots(2,2)

im = density_map(dx_m.ravel(), dx_l.ravel(), ax[0,0])
#im = ax[0,0].scatter(xx, yy, c=z, s=50, edgecolor='')
fig.colorbar(im, ax=ax[0,0], label='Normalized probability')
ax[0,0].set_title(r"$\Delta$ snapshot: " + str(dt), fontsize=24)
ax[0,0].set_xlim([-2,3])
ax[0,0].set_xlabel(r'$\Delta M_{*,f}$', fontsize=24)
ax[0,0].set_ylabel(r'$\Delta \lambda _{R}$', fontsize=24)



from scipy.stats import gaussian_kde

x = dx_m.ravel()
y = dx_l.ravel()

# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy) 
z /= max(z)

z = fg[:,:-dt].ravel()
idx = z.argsort()
xx, yy = x[idx], y[idx]
z = z[idx]
im = ax[0,1].scatter(xx, yy, c=z, s=50, edgecolor='')
fig.colorbar(im, ax=ax[0,1], label='ID')
ax[0,1].set_title(r"$\Delta$ snapshot: " + str(dt), fontsize=24)
ax[0,1].set_xlim([-2,3])
ax[0,1].set_xlabel(r'$\Delta M_{*,f}$', fontsize=24)
ax[0,1].set_ylabel(r'$\Delta \lambda _{R}$', fontsize=24)


z = arrnouts.ravel()
idx = z.argsort()
xx, yy = x[idx], y[idx]
z = z[idx]
im = ax[1,0].scatter(xx, yy, c=z, s=50, edgecolor='')
ax[1,0].set_title(r"$\Delta$ snapshot: " + str(dt), fontsize=24)
ax[1,0].set_xlim([-2,3])
ax[1,0].set_xlabel(r'$\Delta M_{*,f}$', fontsize=24)
ax[1,0].set_ylabel(r'$\Delta \lambda _{R}$', fontsize=24)
fig.colorbar(im, ax=ax[1,0], label='nout')

plt.show()
#%%


