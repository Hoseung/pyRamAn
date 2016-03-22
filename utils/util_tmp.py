# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 17:13:06 2016

@author: hoseung
"""

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


def sigma_clip_ind(c, high, low):
    """
        returns indices of sigma-clipping-safe elements.
    """
    import numpy as np
    ind = (np.mean(c) - np.std(c)*low < c) * (c < np.mean(c) + np.std(c)*high)
    return ind


def mask_outlier(y, low=1.5, high=1.5):
    """
        maks outlier assuming monotonic trend.
    """
    x = np.arange(len(y))

    # linear fitting .. more desirably, a very strong smoothing scheme that can reconstrcut mild curve.
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)

    # extract linear fit
    yy = y - (slope * x + intercept)

    # sigma clipped value = mean of the rest 
    i_good = sigma_clip_ind(yy, low, high)
    yy[~i_good] = np.mean(yy[i_good])

    # add linear fit again
    return yy + (slope * x + intercept)


def smooth(x, beta=5, window_len=20, monotonic=False):
    """ 
    kaiser window smoothing 
    beta = 5 : Similar to a Hamming
    """
    
    if monotonic:
        """
        if there is an overall slope, smoothing may result in offset.
        compensate for that. 
        """
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y=np.arange(len(x)))
        xx = np.arange(len(x)) * slope + intercept
        x = x - xx
    
    # extending the data at beginning and at the end
    # to apply the window at the borders
    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    w = np.kaiser(window_len,beta)
    y = np.convolve(w/w.sum(), s, mode='valid')
    if monotonic: 
         return y[int(window_len)/2:len(y)-int(window_len/2) + 1] + xx#[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    else:
        return y[int(window_len)/2:len(y)-int(window_len/2) + 1]
        #return y[5:len(y)-5]
        
        
def fixed_ind_Lr(gal):
    nnouts = len(gal.nouts)
    ind_reff_fix = np.zeros(nnouts, dtype='i4')

    #print(gal.data['rgal'])
    smooth_r = smooth(mask_outlier(gal.data['rgal'], 1.5, 1.5), 50, monotonic=False)

    # fixed Reff array
    for i in range(nnouts):
        # 1Reff = 5 points
        reff = gal.data['rgal'][i]
        reff_real = smooth_r[i]
        
        try:
            ind_reff_fix[i] = np.round(reff_real/reff * 5) -1
        except:
            pass
    return ind_reff_fix


def smoothed_reff(cat, nout_merger):
    """
    returns "representative" lambda at each nout by assuming monotonic change in Reff. 
    During merger, Reff can fluctuate, and if has no physical meaning to infer Labda at Reff during merger stage. 
    So Reff' is derived by linear interpolating Reffs before and after the merger. 
    
    cat is one galaxy catalog over time.
    """
    import utils.match as mtc
    i_merger = np.where(cat['nout'] == nout_merger)[0]
    ind_lower = 20
    ind_upper = 20
    
    reffs = cat['rgal']
    # left and right values chosen by sigma-clipping
    r_lefts, b, c = scipy.stats.sigmaclip(reffs[max([0,i_merger-ind_lower]):i_merger], sig_lower, sig_upper)
    #print(r_lefts)
    r_left = r_lefts[-1]
    i_left = np.where(reffs == r_left)[0]
    

    r_rights, b,c = scipy.stats.sigmaclip(reffs[i_merger:min([i_merger+ind_upper,len(reffs)])], sig_lower, sig_upper)
    r_right = r_rights[0]
    i_right = np.where(reffs == r_right)[0]

    r_prime = reffs
    #print("chekc")
    #print(r_prime)
    r_prime[i_left : i_right + 1] = np.linspace(r_left, r_right, i_right - i_left + 1)
    return r_prime    

