
# coding: utf-8

# # Measure delta Lambda

# NOTE: Lambda fluctuates, and it fluctuates more as two galaxies get closer.
# It is hard to separate 'normal' stage and 'merging' stage of lambda.
# Measuring L at normal stage may require some fitting algorithm. 

# In[1]:


from utils.util_tmp import *
        
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
            Breaks if i_first == 0.
        """
        # end of galaxy tree = last non-zero position.
        # Note that 'id' can be 0 if phantom. But phantom is a valid datapoint
        i_first_nout = max(np.where(self.data['idx'] > 0)[0])
        print('i_first', i_first_nout)
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
    


# 중간에 비는 것은 아마도 phantom이라 'Orig_halo_id'는 없는 경우일 듯. 
# 그렇기 때문에..! cat을 만들때 idx를 넣어야함!  : 넣었음. (근데 final_ID는 뺐나...?)

# ##### fix Reff of a galaxy

# In[3]:

def func(cluster, ind_lower, ind_upper, sig_lower, sig_upper, nout_fi, bad, fix_ind, mr, dl):
    wdir = '/home/hoseung/Work/data/' + cluster + '/' #'05427/'
    cdir = 'catalog_GM/'

    # Serialize catalogs. -> Only main galaxies

    # main galaxy list
    ## 
    alltrees = ctu.load_tree(wdir, is_gal=True)
    ad = alltrees.data
    tn = ad[ad['nout'] == nout_fi]

    cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout_fi) + '.pickle', 'rb'))
    #idx_all = [tn['id'][tn['Orig_halo_id'] == id_final][0] for id_final in cat['id']]
    last_idx_all = cat['idx'] # idx of galaxies at the last snapshot.
    
    # list of Main progenitor class.
    # feed in alltree.data and idx to initialize MainPrg object.
    mpgs = [MainPrg(ad, idx) for idx in last_idx_all] 

    for nout in range(nout_ini, nout_fi + 1):
        cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout) + '.pickle', 'rb'))
        for mpg in mpgs:
            mpg.set_data(cat, nout) # fill up MainPrg data
        print(nout)

    # load merger galaxy list (geneated by scripts/notebooks/halo/Merter_no_cat.ipynb)
    with open(wdir + 'merger_list.txt', 'rb') as f:
        mgl = np.genfromtxt(f, dtype=[('idx','i8'),('mr','f8'),('nout','i4')])

    mids = mgl['idx']
    mrs = mgl['mr'] 
    nout_mergers = mgl['nout'] 
    i_mergers = nout_fi - nout_mergers 

    gal_idxs = [gal.idxs[0] for gal in mpgs]

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
    #        print("Length", len(gal.nouts))

        if fix_ind:
            ind_reff_fix = fixed_ind_Lr(gal)
            lam = np.zeros(len(ind_reff_fix))

            ind_max = len(gal.data['lambda_arr'][0]) - 1

            for inout, ind in enumerate(ind_reff_fix):#[ind_reff_fix > 0]):
                if ind == 0 : print(ind)
                lam[inout] = gal.data['lambda_arr'][inout][min([ind_max,ind])] # fixed value
        else:
            lam = gal.data['lambda_r']

        x_al = range(max([0,i_merger-ind_lower]), i_merger)
        x_ar = range(i_merger,min([i_merger+ind_upper, len(lam)]))

        al, b1, c1 = scipy.stats.sigmaclip(lam[x_al], sig_lower, sig_upper)
        ar, b2, c2 = scipy.stats.sigmaclip(lam[x_ar], sig_lower, sig_upper)

        if (len(al) > 1) & (len(ar) > 0):
            dl.append(np.median(ar) - np.median(al))
            mr.append(merger_ratio)


# In[4]:

import numpy as np
import scipy.stats
import tree.ctutils as ctu
import matplotlib.pyplot as plt

# Read a single galaxy evolution catalog.
import pickle

clusters = ['36415', '05427', '36413', '39990']
# parameters used for lambda_arr clipping.
ind_upper = 20
ind_lower = 20
sig_upper = 2.0
sig_lower = 2.0

nout_fi = 187

bad = 0

fix_ind=False

mr = []
dl = []

verbose=True

nout_ini = 68
nnouts = nout_fi - nout_ini + 1

cluster = "36415"

wdir = '/home/hoseung/Work/data/' + cluster + '/' #'05427/'
cdir = 'catalog_GM/'

# Serialize catalogs. -> Only main galaxies

# main galaxy list
## 
alltrees = ctu.load_tree(wdir, is_gal=True)
ad = alltrees.data
tn = ad[ad['nout'] == nout_fi]

cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout_fi) + '.pickle', 'rb'))
#idx_all = [tn['id'][tn['Orig_halo_id'] == id_final][0] for id_final in cat['id']]
last_idx_all = cat['idx'] # idx of galaxies at the last snapshot.

# list of Main progenitor class.
# feed in alltree.data and idx to initialize MainPrg object.
mpgs = [MainPrg(ad, idx) for idx in last_idx_all] 


# In[15]:

# Compile catalogs 
for nout in range(nout_ini, nout_fi + 1):
    cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout) + '.pickle', 'rb'))
    for mpg in mpgs:
        mpg.set_data(cat, nout) # fill up MainPrg data
    print(nout, end='\r')

# load merger galaxy list (geneated by scripts/notebooks/halo/Merter_no_cat.ipynb)
with open(wdir + 'merger_list.txt', 'rb') as f:
    mgl = np.genfromtxt(f, dtype=[('idx','i8'),('mr','f8'),('nout','i4')])

mids = mgl['idx']
mrs = mgl['mr'] 
nout_mergers = mgl['nout'] 
i_mergers = nout_fi - nout_mergers 

gal_idxs = [gal.idxs[0] for gal in mpgs]

gal = mpgs[0]
gal.data.dtype

plt.plot(gal.data['xc'])
plt.show()


# In[22]:
fig, axs = plt.subplots(5,5)
axs = axs.ravel()
i = 0

for igal, mid in enumerate(mids):
#gal = mpgs[3]
    if mid not in gal_idxs:
        print("Merger gal {} is not in catalog, skipping".format(mid))
        continue
    else:
        gal = mpgs[np.where(gal_idxs == mid)[0]]

    if len(gal.nouts) < 20:
        continue

    axs[i].plot(gal.data['xc'], 'r')
    gal.clip_non_detection()
    try:
        gal.fill_missing_data()
    except:
        bad = bad + 1
        pass
    axs[i].plot(gal.data['xc'], 'b')
    i +=1
    print(gal.idxs[0])
    if i == 25:
        break


    if verbose: print("Galaxy ID at the final nout, idx = {}, id = {}".format(gal.idxs[0], gal.ids[0]))
    i_merger = i_mergers[igal]  #i_merger = 187 - mgl[igal][2]
    merger_ratio = mrs[igal]

    if verbose: print("nnouts: {}, i_merger {}".format(len(gal.nouts), i_merger))

    if i_merger > len(gal.nouts):
        print("Too short evolution history, aborting..")
        continue

    # fixed Lambda array based on Reff_fix.
#        print("Length", len(gal.nouts))

    if fix_ind:
        ind_reff_fix = fixed_ind_Lr(gal)
        lam = np.zeros(len(ind_reff_fix))

        ind_max = len(gal.data['lambda_arr'][0]) - 1

        for inout, ind in enumerate(ind_reff_fix):#[ind_reff_fix > 0]):
            if ind == 0 : print(ind)
            lam[inout] = gal.data['lambda_arr'][inout][min([ind_max,ind])] # fixed value
    else:
        lam = gal.data['lambda_r']

    x_al = range(max([0,i_merger-ind_lower]), i_merger)
    x_ar = range(i_merger,min([i_merger+ind_upper, len(lam)]))

    al, b1, c1 = scipy.stats.sigmaclip(lam[x_al], sig_lower, sig_upper)
    ar, b2, c2 = scipy.stats.sigmaclip(lam[x_ar], sig_lower, sig_upper)

    if (len(al) > 1) & (len(ar) > 0):
        dl.append(np.median(ar) - np.median(al))
        mr.append(merger_ratio)



# In[8]:

for cluster in clusters:
    func(cluster, ind_lower, ind_upper, sig_lower, sig_upper, nout_fi, bad, fix_ind, mr, dl)

#        except:
#            bad = bad +1
#            pass

print('BAD', bad)


# In[ ]:

# load merger galaxy list (geneated by scripts/notebooks/halo/Merter_no_cat.ipynb)
with open(wdir + 'merger_list.txt', 'rb') as f:
    mgl = np.genfromtxt(f, dtype=[('idx','i8'),('mr','f8'),('nout','i4')])

mids = mgl['idx']
mrs = mgl['mr'] #np.random.random(len(mpgs)) # just as a test.
nout_mergers = mgl['nout'] #np.round(mrs * 50).astype("int") + 137
i_mergers = nout_fi - nout_mergers 

gal_idxs = [gal.idxs[0] for gal in mpgs]

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
#        print("Length", len(gal.nouts))

    if fix_ind:
        ind_reff_fix = fixed_ind_Lr(gal)
        lam = np.zeros(len(ind_reff_fix))

        ind_max = len(gal.data['lambda_arr'][0]) - 1

        for inout, ind in enumerate(ind_reff_fix):#[ind_reff_fix > 0]):
            if ind == 0 : print(ind)
            lam[inout] = gal.data['lambda_arr'][inout][min([ind_max,ind])] # fixed value
    else:
        lam = gal.data['lambda_r']

    x_al = range(max([0,i_merger-ind_lower]), i_merger)
    x_ar = range(i_merger,min([i_merger+ind_upper, len(lam)]))

    al, b1, c1 = scipy.stats.sigmaclip(lam[x_al], sig_lower, sig_upper)
    ar, b2, c2 = scipy.stats.sigmaclip(lam[x_ar], sig_lower, sig_upper)

    if (len(al) > 1) & (len(ar) > 0):
        dl.append(np.median(ar) - np.median(al))
        mr.append(merger_ratio)


# In[12]:

import numpy as np
with open('/home/hoseung/Work/data/01605/m_l.txt', 'rb') as f:
    mrdl  = np.genfromtxt(f)#, dtype=[('idx','i8'),('mr','f8'),('nout','i4')])
mr_1 = mrdl[:,0]
mr_1 = mr_1 + np.random.random(len(mr_1)) * 0.02
dl_1 = mrdl[:,1]
dl_1 = dl_1 + np.random.random(len(dl_1)) * 0.01 - 0.005


# In[ ]:

[mr.append(i) for i in mr_1]
[dl.append(i) for i in dl_1]


# In[28]:

fig, ax = plt.subplots(1)
ax.scatter(mr, dl)
ax.set_ylim([-0.5, +0.5])
ax.set_xlim([0, 10])
ax.set_title("d$\lambda$ vs Merger Mass Ratio")
ax.set_ylabel("d$\lambda$")
ax.set_xlabel("Merger Mass Ratio")
plt.show()


# lambda_mp_gal에서 은하에 대해 np = max 를 i_center로 잡았는데, halo에서 np = max와 다름.....
# 어째 생긴 은하단인가. 허허. 

# I have regularized galaxy evolution data.
# 
# I want to measure rotation parameter before and after a merger.
# I take 20 lambda values before and after the merger, sigma clip outliers, and take median value.
# So dLambda = media(Lambda_after) - median(Lambda_before).
# 
# However, Lambda measurement at 1Reff requires robust Reff measurement, which is very tough during merger events.
# So I smooth Reff evolution history to guess more reasonable Reff values at all points. Following is the procedure.

# In[24]:

len(mr)


# In[ ]:




# In[ ]:




# In[ ]:




# ind_reff_fix points to the Lambda_arr element closest to the fixed Reff at every nout.

# 왜 lambda_r 이랑 lambda_arr[4]랑 다르지? -> 0.5Reff에서 측정했었음..!
# 그림에는 0.5Reff이지만 나머지는 1.0Reff로 쓸래.. 
# 나중에 그림도 1.0으로 바꾸지 뭐.. (lambda_single.py로)

# In[9]:

fig, ax = plt.subplots(1)
ax.plot(lam, 'b-')
ax.plot(lam_fix, 'g--')
plt.show()


# In[11]:

y = gal.data['rgal']
fig, ax = plt.subplots(1)
ax.plot(y)
ax.plot(smooth_r, 'r--')
#ax.plot(smoothed_reff(gal.data['reff'], 102))
ax.plot(gal.data['lambda_r'] * 20)
ax.set_title("Rgal, Rgal_smoothed, and Lambda_r")
plt.show()


# In[19]:

# plot all properties of a galaxy

fig, axs = plt.subplots(4,4)
axs = axs.flatten()
for i, field in enumerate(gal.data.dtype.names[:16]):
    if field == "lambda_arr":
        continue
    axs[i].plot(gal.data[field])
    axs[i].set_ylabel(field)
    
plt.tight_layout()
plt.show()


# In[29]:

# i_merger가 정확해야하는데... 
# Tree에서 주는 merger는 final coalescence일 가능성이 높음. 
print(np.mean(al), np.mean(ar))
print(np.median(al), np.median(ar))


fig, ax = plt.subplots(1)
ax.plot(x_al, lam[x_al], 'r')
ax.plot(x_ar, lam[x_ar], 'b')
ax.axhline(np.mean(al), color='r')
ax.axhline(np.mean(ar), color='b')
#ax.plot(smoothed_reff(gal.data['reff'], 102))
ax.set_title("Rgal, Rgal_smoothed, and Lambda_r")
plt.show()


# I can measure dL in this way, but am I following the right galaxy? is the tree right?
