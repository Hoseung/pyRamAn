# coding: utf-8
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

def density_map(x, y, sort=True):
    from scipy.stats import gaussian_kde
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy) 
    z /= max(z)

    idx = z.argsort()    
    xx, yy = x[idx], y[idx]
    z = z[idx]
    
    return xx, yy, z


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

        
class MainPrg():
    import tree.ctutils as ctu
    import numpy as np
    
    def __init__(self, treedata, final_gal, nout_ini=None, nout_fi=None):

        temp_tree = ctu.extract_main_tree(treedata, final_gal)
        if nout_ini == None:
            nout_ini = min(temp_tree['nout'])
        if nout_fi == None:
            nout_fi = max(temp_tree['nout'])            
            
        self.nouts = np.arange(nout_fi, nout_ini -1, -1)
        self.idxs = temp_tree['id'] # nout_ini, nout_fi consideration needed.
        self.ids = temp_tree['Orig_halo_id']
        self.data = None
    
    def set_data(self, cat, nout):
        """
        compile data from catalogs.
        """
        if nout in self.nouts:
            # Declare self.data first if there isn't.
            if self.data == None:
                self.data = np.zeros(len(self.nouts), dtype=cat.dtype)
            inow = self.nouts == nout
            a = np.where(cat['idx'] == self.idxs[inow])[0]
            if len(a) > 0:
                self.data[inow] = cat[a]        
            else:
                pass
                #print(self.ids[inow],cat['id'])
        else:
            pass
            #print("No {} in the catalog".format(nout))
            
    def clip_non_detection(self):
        # end of galaxy tree = last non-zero position.
        # Note that 'id' can be 0 if phantom. But phantom is a valid datapoint
        i_first_nout = max(np.where(self.data['idx'] > 0)[0])
        #print('i_first', i_first_nout)
        # then, only [0-i_first_nout] are valid.
        # earlier then 187 - 91-th are zero. so get rid of them.
        self.data = self.data[:i_first_nout].copy()
        self.nouts = self.nouts[:i_first_nout].copy()
        self.ids = self.ids[:i_first_nout].copy()
        self.idxs = self.idxs[:i_first_nout].copy()
        
    def fill_missing_data(self):
        assert (self.ids[-1] != 0)
        
        # loop over all fields except id, index, and non-physical entries.
        i_bad = np.where(self.data['idx'] == 0)[0]
        for field in self.data.dtype.names:
            # do not modify index and id fields.
            if field in ["index", "id", "idx"]:
                continue
            arr = self.data[field] # it's a view.

            for i_b in i_bad:
                # neighbouring array might also be empty. Search for closest valid element.
                # left point
                i_l = i_b - 1
                while(i_l in i_bad):
                    i_l = i_l - 1
                # right point
                i_r = i_b + 1
                while(i_r in i_bad):
                    i_r = i_r + 1

                arr[i_b] = (arr[i_b -1] + arr[i_b +1])/2.
    


# In[2]:

def fixed_ind_Lr(gal):
    nnouts = len(gal.nouts)
    ind_reff_fix = np.zeros(nnouts, dtype='i4')

    #print(gal.data['rgal'])
    smooth_r = smooth(mask_outlier(gal.data['rgal'], 1.5, 1.5), 50, monotonic=False)

    # fixed Reff array
    for i in range(nnouts):
        # 1Reff = 5 points
        reff_real = smooth_r[i]
        reff = gal.data['rgal'][i]
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


# In[3]:

import numpy as np
import scipy.stats
import tree.ctutils as ctu
import matplotlib.pyplot as plt

# Read a single galaxy evolution catalog.
import pickle


# In[4]:

clusters = ['10002', '04466', '17891', '36415', '35663', '06098', '07206',\
            '49096', '39990', '36413', '01605', '05427'][:]
# parameters used for lambda_arr clipping.
ind_upper = 20
ind_lower = 20
sig_upper = 2.0
sig_lower = 2.0

nout_ini = 70
nout_fi = 187

bad = 0


# In[ ]:

base = '/data1/good/'
cdir = ['catalog/', 'easy/', 'catalog_GM/'][1]


verbose=True

ngals_tot = 0

for cluster in clusters:
    wdir = base + cluster + '/'
    # main galaxy list
    cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout_fi) + '.pickle', 'rb'))
    ngals_tot = ngals_tot + len(cat['idx'])


nnouts = nout_fi - nout_ini + 1

mpgs = []
for cluster in clusters:
    print(cluster)
    wdir = base + cluster + '/'

    # Serialize catalogs. -> Only main galaxies
    # main galaxy list
    alltrees = ctu.load_tree(wdir, is_gal=True)
    ad = alltrees.data
    tn = ad[ad['nout'] == nout_fi]

    cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout_fi) + '.pickle', 'rb'))
    #idx_all = [tn['id'][tn['Orig_halo_id'] == id_final][0] for id_final in cat['id']]
    idx_all = cat['idx']

    mpg_tmp = []
    for i, idx in enumerate(idx_all):
        mpg_tmp.append(MainPrg(ad, idx))
        print(i, idx)
#    mpg_tmp =[MainPrg(ad, idx) for idx in idx_all]
    for nout in range(nout_ini, nout_fi + 1):
        cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout) + '.pickle', 'rb'))
        for gal in mpg_tmp:
            gal.set_data(cat, nout)
        print(nout)

    while len(mpg_tmp) > 0:
        mpgs.append(mpg_tmp.pop())

    
with open('main_prgs_GM.pickle', 'wb') as f:
    pickle.dump(mpgs, f)

