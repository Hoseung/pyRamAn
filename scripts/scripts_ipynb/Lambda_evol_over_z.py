
# coding: utf-8

# In[1]:

import pickle
import matplotlib.pyplot as plt
import numpy as np
from analysis.evol_lambda import MainPrg


# Gaussian KDE 
# 
# Look at this article for the idea
# http://www.mglerner.com/blog/?p=28
# and this for comparison between implementations.
# https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
# 

# In[2]:

from analysis.all_plot_modules import *
import analysis.Major_Minor_accretion as mma


# In[3]:

wdir = '/home/hoseung/Work/data/'
raw_mpgs = False
save=True
if raw_mpgs:
    #mpgs = pickle.load(open(wdir + "all_prgs/main_prgs_ALL.pickle", "rb"))
    mpgs = pickle.load(open(wdir + "all_prgs/main_prgs_final_augmented_5_10_0.5_0.5_1.0_37_0.01_filtered_.pickle", "rb"))

    for gal in mpgs:
        #print(sum(gal.data["reff"] > 0),len(smooth(gal.data["reff"])))
        #smoothed_reff = smooth(gal.data["reff"])
        gal.fill_missing_data()
        #print(1,min(gal.nouts))
        gal.clip_non_detection()
        #print(2,min(gal.nouts))
        gal.smoothed_lambda_org = mma.smooth(gal.data["lambda_r"], window_len=15)[:-1]
        #print(3,min(gal.nouts))
        gal.smoothed_r = mma.smooth(gal.data["reff"], window_len=15)[:-1]
        #print(4,min(gal.nouts))
        gal.smoothed_lambda = mma.smooth(l_at_smoothed_r(gal, npix_per_reff=5), window_len=15)[:-1]
        #print(5,min(gal.nouts))
        if save:
            pickle.dump(mpgs, open(wdir + "main_prgs_final_augmented_5_10_0.5_0.5_0.5_37_0.01_filtered.pickle", "wb"))
else:
    #mpgs = pickle.load(open("main_prgs_final_augmented_5_10_0.5_0.5_0.5_37.pickle", "rb"))
    mpgs = pickle.load(open(wdir + "main_prgs_final_augmented_5_10_0.5_0.5_0.5_37_0.01_filtered_.pickle", "rb"))


# In[5]:

def den_lambda_evol(mpgs, nout_ini, nout_fi,
                     wdir_info='./',
                     density = "none",
                     sizeOfFont=9):


    nnouts = nout_fi - nout_ini + 1
    ngals_tot = len(mpgs)
    lambda_evol_all = np.zeros([ngals_tot, nnouts])

    # Starting nouts of galaxies are different. 
    # But all end at nout = 187.
    for igal, gal in enumerate(mpgs):
        for inout, nout in enumerate(gal.nouts):
            lambda_evol_all[igal][nout - nout_ini] = gal.data['lambda_r'][inout]


    # change ticks
    zreds=[]
    aexps=[]

    # Only nout_ini < nout < nout_fi values are taken.
    # So are the x ticks.
    import load
    for nout in range(nout_ini, nout_fi + 1):
        info = load.info.Info(nout=nout, base=wdir_info, load=True)
        aexps.append(info.aexp)
        zreds.append(info.zred)
    aexps = np.array(aexps)
    zreds = np.array(zreds)

#    modify_ticks1(zreds, aexps, axs[2], nout_ini, nout_fi, fontsize=fontsize_tick_label)

    xx = np.tile(np.arange(nnouts), ngals_tot)
    all_data = lambda_evol_all.ravel()

    ind_ok = all_data > 0.01
    lambda_range=[0.01, 0.8]
    xx,yy,z = density_map(xx[ind_ok], all_data[ind_ok])
    
    return xx,yy,z, zreds, aexps


def plot_lambda_evol(xx,yy,z,axs,cmap ="jet", img_scale=1.5):
    fontsize_ticks = 6 * img_scale
    fontsize_tick_label = 8 * img_scale
    fontsize_legend = 5 * img_scale

    for ax in axs:
        ax.scatter(xx, yy, c=z, s=50, edgecolor='', cmap=cmap, rasterized=True)
#            yticks_ok=[0.0, 0.2, 0.4, 0.6, 0.8]
#            ax.set_ylim([-0.05, 0.9])
#            ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
#            ax.set_yticklabels([str(yy) for yy in yticks_ok])
#            ax.set_ylabel(r"$\lambda_{R_{eff}}$", fontsize=fontsize_tick_label, family="Liberation Sans")
#            ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)

#    axs[2].tick_params(axis='x', which='major', labelsize=fontsize_ticks)
#    return zreds, aexps


# In[7]:

nout_ini = 37
fig, axs = plt.subplots(3, sharex=True)
fig.set_size_inches(4.75, 8)
plt.subplots_adjust(hspace=0.01)


# In[12]:

xx,yy,z, zreds, aexps = den_lambda_evol(mpgs, nout_ini, 187,
             wdir_info = wdir + '29172/',
             density="kernel")

#pickle.dump((xx,yy,z,zreds,aexps), open("lambda_evol_xxyyz.pickle", "wb"))


# In[8]:

xx,yy,z,zreds, aexps = pickle.load(open(wdir+"lambda_evol_xxyyz.pickle", "rb"))

plot_lambda_evol(xx,yy,z, axs, cmap = "Blues")  # 37 : z = 3.0


# In[21]:

#pickle.dump(fig, open("lambda_evol_background.pickle", "wb"))
#fig = pickle.load(open("lambda_evol_background.pickle", "rb"))


# In[10]:

plt.show()


# In[9]:

plt.savefig("figs/fig4_bare.png", dpi=200)
plt.savefig("figs/fig4_bare.pdf")


# In[10]:

len(mpgs)


# In[9]:

selected_cluster = 29176

#fig, ax = plt.subplots()
#selected_gals = [40, 75] # major mergers
#gals = [i for i, gal in enumerate(mpgs) if (gal.cluster == selected_cluster and gal.ids[0] in selected_gals)]

i_major =[13, 34, 39, 43, 44, 51, 116, 219, 220, 728, 742]
#gals = [gal for i, gal in enumerate(mpgs) if i in [340, 342, 349]]
gals = [mpgs[728], mpgs[226], mpgs[742]]
plot_major(gals, axs[0],
           suptitle="Major Mergers",
           img_scale=1.5,
           arrow_scale=20) # arrow_scale = 50 for png, 20 for vectors.

#selected_gals = [17] # minor mergers
#gals = [gal for gal in mpgs if (gal.cluster == selected_cluster and gal.ids[0] in selected_gals)]
gals = []
gals.append(mpgs[285])
gals.append(mpgs[507])
plot_minor(gals, axs[1],
              suptitle="Minor Mergers",
              style="stream",
              img_scale=1.5,
              annotate="(B) ",
              arrow_scale=20)

selected_gals = [45, 50, 98, 130] # no mergers
gals = [gal for gal in mpgs if (gal.cluster == selected_cluster and gal.ids[0] in selected_gals)]
plot_rest(gals, axs[2],
          suptitle="No Mergers",
          style="stream",
          img_scale=1.5,
          annotate="(C) ",
          arrow_scale=20)

plt.tight_layout()


#Fix font

from matplotlib import rc, font_manager
sizeOfFont=9
fontProperties = {'family':'Liberation Sans',
                  'weight' : 'normal', 'size' : sizeOfFont}
ticks_font = font_manager.FontProperties(family='Liberation Sans', style='normal',
               size=sizeOfFont, weight='normal', stretch='normal')
rc('text', usetex=True)
#rc('font',**fontProperties)

fontsize_ticks = 8
fontsize_tick_label = 10
for ax in axs:
    ax.set_yticklabels(ax.get_yticks(), fontProperties)
    yticks_ok=[0.0, 0.2, 0.4, 0.6, 0.8]
    ax.set_ylim([-0.05, 0.9])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_yticklabels([str(yy) for yy in yticks_ok])
    ax.set_ylabel(r"$\lambda_{R_{eff}}$", fontsize=fontsize_tick_label, family="Liberation Sans")
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    
axs[2].tick_params(axis='x', which='major', labelsize=fontsize_ticks)
axs[2].set_xticklabels(axs[2].get_xticks(), fontProperties)

nout_ini = 37
nout_fi = 187
modify_ticks1(zreds, aexps, axs[2], nout_ini, nout_fi, fontsize=9)


fname_base = wdir+"figs/lambdar_evol_test"
plt.savefig(fname_base + "smooth_"+".pdf")
#plt.savefig(fname_base + "smooth_"+".eps")
#plt.savefig(fname_base + "smooth_"+".svg")
plt.savefig(fname_base + "smooth_"+".png", dpi=200)
plt.savefig(fname_base + "smooth_Hires"+".png", dpi=400)


# In[ ]:




# In[ ]:



