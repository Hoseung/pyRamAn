
# coding: utf-8

# # Measure delta Lambda

# NOTE: Lambda fluctuates, and it fluctuates more as two galaxies get closer.
# It is hard to separate 'normal' stage and 'merging' stage of lambda.
# Measuring L at normal stage may require some fitting algorithm. 

# In[1]:

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
        """
            remove empty data points at early universe where galaxy were not 
            massive enough to be detected. 
            Assuming that "end of galaxy tree" = "last non-zero position".
            
            Bug: Breaks if i_first == 0.
            
        """
        # Note that 'id' can be 0 if phantom. But phantom is a valid datapoint
        i_first_nout = max(np.where(self.data['idx'] > 0)[0])
        
        
        
        i_first_nout = max([i_first_nout, 10]) # 임시로... 
        
        
        
        print('i_first', i_first_nout)
        # then, only [0 ~ i_first_nout] are valid.
        # earlier then 187 - 91-th are zero. get rid of them.
        self.data = self.data[:i_first_nout].copy()
        self.nouts = self.nouts[:i_first_nout].copy()
        self.ids = self.ids[:i_first_nout].copy()
        self.idxs = self.idxs[:i_first_nout].copy()
        
    def fill_missing_data(self):
        """
            detect missing data and fill the hole by interpolating adjascent data points.
        """
        
        assert (self.ids[-1] != 0)
        
        i_bad = np.where(self.data['idx'] == 0)[0]
        
        # loop over all physical quantiies (except id, index, and so on).
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
    


# 중간에 비는 것은 아마도 phantom이라 'Orig_halo_id'는 없는 경우일 듯. 
# 그렇기 때문에..! cat을 만들때 idx를 넣어야함!  : 넣었음. (근데 final_ID는 뺐나...?)

# ##### fix Reff of a galaxy

# In[2]:

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



def get_mr_dl(cluster, mr, dl):
   
    wdir = './' + cluster + '/' #'05427/'
    # Serialize catalogs. -> Only main galaxies
    # main galaxy list
    alltrees = ctu.load_tree(wdir, is_gal=True)
    ad = alltrees.data
    tn = ad[ad['nout'] == nout_fi] # tree now

    cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout_fi) + '.pickle', 'rb'))
    last_idx_all = cat['idx'] # idx of galaxies at the last snapshot.

    # list of Main progenitor class.
    #   feed in alltree.data and idx to initialize MainPrg object.
    mpgs = [MainPrg(ad, idx) for idx in last_idx_all] 
    
    # Compile catalogs ##################
    for nout in range(nout_ini, nout_fi + 1):
        cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout) + '.pickle', 'rb'))
        for mpg in mpgs:
            mpg.set_data(cat, nout) # fill up MainPrg data
        print(nout, end='\r')
        
     
    # load merger galaxy list (geneated by scripts/notebooks/halo/Merter_no_cat.ipynb) ##################
    with open(wdir + 'merger_list.txt', 'rb') as f:
        mgl = np.genfromtxt(f, dtype=[('idx','i8'),('mr','f8'),('nout','i4')])

    mids = mgl['idx'] # Merger-galaxy IDs
    mrs = mgl['mr'] # Merger ratios
    nout_mergers = mgl['nout'] # merger epoch
    i_mergers = nout_fi - nout_mergers # merger epoch as nout index.

    gal_idxs = [gal.idxs[0] for gal in mpgs] # idx of all galaxies from cluster catalog.

    
    bad = 0
    
    #1. merger list에 있는 은하들에 대해서 은하 catalog를 업데이트 한다. 
    #2. 은하의 catalog가 충분히 긴지 확인
    #3. Rvir이 요동치는 경우 Reff를 보정. (요동치기 전, 후의 값을 linear interpolate)
    #4. 전-후의 lambda를 측정해서 dl vs mr 측정. 
    for igal, mid in enumerate(mids):
    #gal = mpgs[3]
        if mid not in gal_idxs:
            print("Merger gal {} is not in catalog, skipping".format(mid))
            continue
        else:
            gal = mpgs[np.where(gal_idxs == mid)[0]]

        if len(gal.nouts) < 20:
            continue
        gal.clip_non_detection()
        try:
            gal.fill_missing_data()
        except:
            bad = bad + 1
            pass


        if verbose: print("Galaxy ID at the final nout, idx = {}, id = {}".format(gal.idxs[0], gal.ids[0]))

        i_merger = i_mergers[igal]  #i_merger = 187 - mgl[igal][2]
        merger_ratio = mrs[igal]

        if verbose: print("nnouts: {}, i_merger {}".format(len(gal.nouts), i_merger))

        if i_merger > len(gal.nouts):
            print("Too short evolution history, aborting..")
            continue

        # fixed Lambda array based on Reff_fix.
        if fix_ind:
            ind_reff_fix = fixed_ind_Lr(gal)
            lam = np.zeros(len(ind_reff_fix))

            ind_max = len(gal.data['lambda_arr'][0]) - 1

            for inout, ind in enumerate(ind_reff_fix):#[ind_reff_fix > 0]):
                if ind == 0 : print(ind)
                lam[inout] = gal.data['lambda_arr'][inout][min([ind_max,ind])] # fixed value
        else:
            lam = gal.data['lambda_r']

        x_al = range(max([0,i_merger-ind_lower]), i_merger) # nout before merger
        x_ar = range(i_merger,min([i_merger+ind_upper, len(lam)])) # nout after merger

        # representative value of lambda BEFORE the merger
        al, b1, c1 = scipy.stats.sigmaclip(lam[x_al], sig_lower, sig_upper) 

        # representative value of lambda AFTER the merger
        ar, b2, c2 = scipy.stats.sigmaclip(lam[x_ar], sig_lower, sig_upper)

        if (len(al) > 1) & (len(ar) > 0):
            dl.append(np.median(ar) - np.median(al))
            mr.append(merger_ratio)
        else:
            print("error in measuring lambda")
    

# In[27]:

import numpy as np
import scipy.stats
import tree.ctutils as ctu
import matplotlib.pyplot as plt

# Read a single galaxy evolution catalog.
import pickle


# parameters used for lambda_arr clipping.
ind_upper = 20
ind_lower = 20
sig_upper = 2.0
sig_lower = 2.0

nout_fi = 187


fix_ind=False

verbose=True

nout_ini = 68
nnouts = nout_fi - nout_ini + 1

cdir = 'catalog_GM/'

clusters = ['17891','35663', '39990', '36413', '01605', '36415', '05427']#17891
mr = [] # mass ratio
dl = [] # delta lambda
fig,ax = plt.subplots(1)
for cluster in clusters:
    print(" %%%%%%%%%% " , cluster, " %%%%%%%%%% ")
    get_mr_dl(cluster, mr, dl)
    ax.scatter(mr,dl)

pickle.dump([dl, mr],open('dlmr.pickle', 'wb'))
# In[31]:

#plt.scatter(mr,dl)
plt.savefig("dlmr.png")
