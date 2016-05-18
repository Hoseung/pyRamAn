
# coding: utf-8

# ## Measure contribution of Major merger, minor merger, and smooth accretion for only the 'safe' samples. 
# Because tree bad link more likely occur at major merger events, I guess the 'safe' samples have less major mergers than the total sample. 
# 
# * 29176 catalog is not complete.
# change to 10002

# ## import analysis.Major_Minor_accretion as mma
# 

# In[12]:

import utils.sampling as smp
import matplotlib.pyplot as plt
import tree
import pickle
import tree.halomodule as hmo
import numpy as np
import draw
import load
import analysis.evol_lambda as evl
import analysis.Major_Minor_accretion as mma
import tree.ctutils as ctu
import utils.match as mtc


nout_fi = 187
nout_ini = 57
wdir = './'
cdir = 'easy/'
dir_out = wdir 
suffix = "10002"

info = load.info.Info(nout=187, load=True)


# Largest halo in the refinement region is not the main cluster.
# Is the largest NP halo the main cluster? 
# To check it, color halos in NP.
# 
# No, max_np galaxy/halo is the main galaxy/halo.
# But 'rvir' value is wrong.
# 
# and fixed.

# In[13]:

gg = hmo.Halo(nout=187, base=wdir, load=True, is_gal=True)
hh = hmo.Halo(nout=187, base=wdir, load=True, is_gal=False)

good_gals = mma.close_gals(hh, gg, rscale=3)


# In[9]:

# cluster 29176 ok_gals=[3, 5, 6, 9, 10, 11, 15, 16, 18, 20, 21, 22, 25, 26, 31, 32, 34, 37, 39, 43, 45, 49, 54, 55, 62, 63, 67, 68, 71]
ok_gals =[2, 3,4, 6, 7, 10, 13, 18, 19, 22, 23, 28, 34,
          40, 41, 42, 47, 48, 50, 51, 52, 54, 55, 59,
          65, 66, 68, 73, 75, 81, 82, 83, 88, 89, 101, 102,
          108, 113, 128, 131, 132, 135, 144, 148, 151, 157, 158,
          170, 171, 178, 179,182,208, 214, 225, 421] # 10002
# Galaxies with good main progenitor list.
# But some galaxies have up-and-down mass evolution histories. 
# By looking at merger tree plots, I can roughly tell 
# major merger galaixies, minor merger galaxies, smooth accretion galaxies.


# Load and compile catalogs

# In[6]:

alltrees = ctu.load_tree(wdir, is_gal=True)
# Serialize catalogs. -> Only main galaxies
# main galaxy list

ad = alltrees.data
tn = ad[ad['nout'] == nout_fi]

cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout_fi) + '.pickle', 'rb'))
# Only good galaxies are saved in the catalogue.
# But there is no mass cut in good_gals list.
# And match_list_ind doest not handle cases where the array a is shorter than array b.
# Need another matching funciton that automatically excludes missing values.
good_gals = np.intersect1d(cat['id'], good_gals)
ind_ok = mtc.match_list_ind(cat['id'], good_gals, allow_swap=False)
cat = cat[ind_ok]

# Some galaxies are not in the tree. So can't follow their main progenitors.
cat = cat[~np.isnan(cat['idx'])]
idx_all = cat['idx'].astype(int) # fix so that idx is int from the beginning. 
#print(good_gals)


mpgs = mma.compile_mpgs(alltrees, idx_all, wdir=wdir, cdir=cdir, nout_ini=nout_ini, nout_fi=nout_fi)


# Merger epochs

# In[7]:

mma.find_merger_epochs(alltrees, idx_all, mpgs, nout_ini=nout_ini)


# ## Mass evolution, lambda evolution, smoothed lambda evolution plots

# ### Measure delta lambda

# In[8]:

mma.measure_delta_lambda(mpgs, dt=7, nout_ini=nout_ini)


# In[17]:

from matplotlib.backends.backend_pdf import PdfPages

fig, ax = plt.subplots(2,2, sharex=True)
plt.subplots_adjust(hspace=0.001)

with PdfPages(wdir + 'multipage_pdf.pdf') as pdf:
    for gal in mpgs:
        ind_nout = gal.nouts > nout_ini
        gal.nouts = gal.nouts[ind_nout]
        gal.data = gal.data[ind_nout]
        #smoothed_lambda = medfilt(gal.data['lambda_r'], kernel_size=5)
        
        ax[0,0].scatter(gal.nouts, np.log10(gal.data['mstar']))
        ax[0,0].set_xlim([50,190])
        ax[0,0].set_title(str(gal.ids[0]) + ", " + str(gal.idxs[0]))
        #ax[0].set_ylim([8,13])
        ax[1,0].plot(gal.nouts, gal.smoothed_lambda, 'r-')
        ylim_mid = ax[1,0].get_ylim() # match ylim of original and smoothed plot
        ax[0,1].plot(gal.nouts, gal.data['lambda_r'], 'r-')
        ax[0,1].set_ylim(ylim_mid)
        
        
        
        if gal.merger is not None:
            delta_lambda =[]
            for mr, xx, x2 in zip(gal.merger.mr, gal.merger.nout, gal.merger.nout_ini):
                ax[0,0].axvline(xx, linestyle=':')
                ax[0,0].annotate("{:.1f}".format(mr), xy=(xx,0.8))
                ax[0,1].axvline(xx, linestyle=':')
                ax[0,1].annotate("{:.1f}".format(mr), xy=(xx,0.8))
                ax[1,0].axvline(xx, linestyle=':')
                ax[1,0].axvline(x2, linestyle=':', c='g')
                
                i_nout = np.where(gal.nouts == xx)[0]
                #print(max([0, i_nout - dt]))
                #print(i_nout)
                
                # nout[0] = 187
                #lambda_after = np.average(smoothed_lambda[max([0, i_nout - dt]) : i_nout])
                #lambda_before = np.average(smoothed_lambda[i_nout:min([len(gal.data), i_nout + dt])])
                #delta_lambda.append(lambda_after - lambda_before)
            #gal.merger.delta = np.array(delta_lambda)
        pdf.savefig()
        ax[0,0].clear()
        ax[0,1].clear()
        ax[1,0].clear()

#plt.show()
plt.close()


# ### simple plot

# In[9]:

mma.Maj_min_acc_ratio(mpgs)
pickle.dump(mpgs, dir_out + "mpgs" + suffix + ".pickle")


# In[10]:

fig, ax = plt.subplots(2)
for igal, gal in enumerate(mpgs):
    if gal.merger is not None:
        for mr, xx, delta in zip(gal.merger.mr, gal.merger.nout, gal.merger.delta):
            ax[0].scatter(mr, delta)
    ax[1].scatter([1],gal.dlM) 
    ax[1].scatter([3],gal.dlm) 
    ax[1].scatter([5],gal.dlo)
#plt.show()
plt.savefig(dir_out + "mma_simple" + suffix + ".png")


# ### Draw a vilon plot of Major/Minor/accretion contributio ratio

# In[11]:

dlt_all=[]
dlo_all=[]
dlM_all=[]
dlm_all=[]
for igal, gal in enumerate(mpgs):
    dlt_all.append(gal.dlt)
    dlo_all.append(gal.dlo)
    dlM_all.append(gal.dlM)
    dlm_all.append(gal.dlm)
    

data = [dlM_all, dlm_all, dlt_all]
pos = [1,2,3]
fig, ax = plt.subplots()
ax.violinplot(data, pos, points=20, widths=0.3,
                      showmeans=True, showextrema=True, showmedians=True)
ax.set_ylim([-0.8, 0.3])
ax.tick_params(labelsize=18)
plt.savefig(dir_out + "mma_violin" + suffix + ".png")
#plt.show()


# In[ ]:



