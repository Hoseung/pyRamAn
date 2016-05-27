
# coding: utf-8

# In[1]:

def multipage_plots(mpgs, nout_ini=37, wdir='./', suffix=""):
    from matplotlib.backends.backend_pdf import PdfPages

    fig, ax = plt.subplots(2,2, sharex=True)
    plt.subplots_adjust(hspace=0.001)

    with PdfPages(wdir + 'MMA_plots' + suffix +'.pdf') as pdf:
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


# In[2]:

def plot_violin(mpgs, mstar_limit=1e10, suffix="", wdir='./'):
    dlt_all=[]
    dlo_all=[]
    dlM_all=[]
    dlm_all=[]
    for igal, gal in enumerate(mpgs):
        # Only large galaxy at the final snapshot
        if gal.data['mstar'][0] > mstar_limit:
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


# In[3]:

def plot_simple(mpgs, wdir='./', suffix=""):
    fig, ax = plt.subplots(2)
    for igal, gal in enumerate(mpgs):
        if gal.merger is not None:
            for mr, xx, delta in zip(gal.merger.mr, gal.merger.nout, gal.merger.delta):
                ax[0].scatter(mr, delta)
        ax[1].scatter([1],gal.dlM) 
        ax[1].scatter([3],gal.dlm) 
        ax[1].scatter([5],gal.dlo)
    #plt.show()
    plt.savefig(wdir + "mma_simple" + suffix + ".png")


# In[18]:

def plot_hist(mpgs, wdir=wdir, suffix=""):
    all_mr = []
    all_delta =[]
    for gal in mpgs:
        if gal.merger is not None:
            for mr, delta in zip(gal.merger.mr, gal.merger.delta):
                all_mr.append(mr)
                all_delta.append(delta)

    all_mr = np.array(all_mr)
    all_delta = np.array(all_delta)

    plt.hist(all_delta[all_mr < 3], alpha=0.3, range=[-0.7, 0.5],
             facecolor='none', edgecolor='r', linestyle='dashed',
             lw=3)
    plt.hist(all_delta[(all_mr < 10) * (all_mr > 3)], alpha=0.3, range=[-0.7, 0.5],
             facecolor='none', edgecolor='g', linestyle='solid', lw=3)
    plt.hist(all_delta[all_mr > 10], alpha=0.3, range=[-0.7, 0.5],
             facecolor='none', edgecolor='b', linestyle="-.", lw=3)
    plt.savefig(wdir + 'MMA_hist' + suffix +'.png')
    
    print(len(all_mr))


# In[8]:

def filter_small_mergers(mm, window=7):
    if mm is None:
        return
    
    if len(mm.nout) == 1:
        return

    
    neighbor=np.zeros(len(mm.mr))
    # Tag mergers with neighbor
    delimeter = 1
    for i in range(len(mm.mr)-1):
        #print(mm.nout[i], mm.nout[i+1], delimeter)
        neighbor[i]=delimeter
        if abs(mm.nout[i] - mm.nout[i+1]) > window:
            delimeter +=1

        
        #print(neighbor)
    # last element
    neighbor[-1] = delimeter
    #if abs(mm.nout[-2] - mm.nout[-1]) < window:
    #    neighbor[-1] = delimeter
    #else:
    #    neighbor[-1] = delimeter + 1

    #print(mm.nout)
    #print(neighbor)

    ind_ok = []
    # neighbor may have 0,1, or may start from 2.
    for imerger in range(delimeter+1):
        ii = np.where(neighbor == imerger)[0]
        if len(ii) > 0:
            ind_ok.append(ii[np.argmin(mm.mr[ii])])

    #print("finally", neighbor[ind_ok],mm.mr, mm.mr[ind_ok])
    mm.nout = mm.nout[ind_ok]
    mm.mr = mm.mr[ind_ok]
    mm.nout_ini = mm.nout_ini[ind_ok]


# In[11]:

def multipage_plots(mpgs, nout_ini=37, wdir='./', suffix="", dt_after = 10, dt_before = 7):
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.signal import medfilt

    fig, ax = plt.subplots(1,2, sharex=True)
    plt.subplots_adjust(hspace=0.001)
    fig.set_size_inches(8,4)

    with PdfPages(wdir + 'MMA_plots' + suffix +'.pdf') as pdf:
        for gal in mpgs:
            ind_nout = gal.nouts > nout_ini
            gal.nouts = gal.nouts[ind_nout]
            gal.data = gal.data[ind_nout]

            ax[0].scatter(gal.nouts, np.log10(gal.data['mstar']))
            ax[0].set_xlim([50,190])
            ax[0].set_title(str(gal.ids[0]) + ", " + str(gal.idxs[0]))
            ax[1].plot(gal.nouts, gal.smoothed_lambda, 'black')

            if gal.merger is not None:
                delta_lambda =[]

                ind_nout = gal.nouts > nout_ini
                gal.nouts = gal.nouts[ind_nout]
                gal.data = gal.data[ind_nout]
                gal.smoothed_lambda = medfilt(gal.data['lambda_r'], kernel_size=5)

                for mr, xx, x2 in zip(gal.merger.mr, gal.merger.nout, gal.merger.nout_ini):
                    i_nout = np.where(gal.nouts == xx)[0]
                    iini_nout = np.where(gal.nouts == x2)[0]

                    nouts_after = gal.nouts[max([0, i_nout - dt_after]) : i_nout]
                    l_after = gal.smoothed_lambda[max([0, i_nout - dt_after]) : i_nout]
                    lambda_after = np.average(l_after)
                    
                    ax[1].plot(nouts_after,l_after, 'g-')
                    nn = range(min(nouts_after) - 5, max(nouts_after) + 5)
                    ax[1].plot(nn, [lambda_after]*len(nn), "g:")

                    l_before = gal.smoothed_lambda[iini_nout:min([len(gal.data), iini_nout + dt_before])]
                    nouts_before = gal.nouts[iini_nout:min([len(gal.data), iini_nout + dt_before])]
                    lambda_before = np.average(l_before)
                    ax[1].plot(nouts_before,l_before, 'r-')
                    nn = range(min(nouts_before) - 5, max(nouts_before) + 5)
                    ax[1].plot(nn, [lambda_before]*len(nn), "r:")

                    delta_lambda.append(lambda_after - lambda_before)
                gal.merger.delta = np.array(delta_lambda)

                for mr, xx, x2 in zip(gal.merger.mr, gal.merger.nout, gal.merger.nout_ini):
                    ax[0].axvline(xx, linestyle=':')
                    ax[0].annotate("{:.1f}".format(mr), xy=(xx,0.8))
                    ax[1].axvline(xx, linestyle=':')
                    ax[1].annotate("{:.1f}".format(mr), xy=(xx,0.8))
                    ax[1].axvline(xx, linestyle=':')
                    ax[1].axvline(x2, linestyle=':', c='g')

                    i_nout = np.where(gal.nouts == xx)[0]
            
            pdf.savefig()
            ax[0].clear()
            ax[1].clear()

    #plt.show()
    plt.close()


# ## Measure contribution of Major merger, minor merger, and smooth accretion for only the 'safe' samples. 
# Because tree bad link more likely occur at major merger events, I guess the 'safe' samples have less major mergers than the total sample. 
# 
# * 29176 catalog is not complete.
# change to 10002

# ## import analysis.Major_Minor_accretion as mma
# 

# In[4]:

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
cdir = 'easy_new/'
suffix = ["10002", "04466", "29172"][2]
wdir = './' + suffix + '/'
dir_out = wdir 
# "10002"

info = load.info.Info(base= wdir, nout=187, load=True)


# Largest halo in the refinement region is not the main cluster.
# Is the largest NP halo the main cluster? 
# To check it, color halos in NP.
# 
# No, max_np galaxy/halo is the main galaxy/halo.
# But 'rvir' value is wrong.
# 
# and fixed.

# In[5]:

gg = hmo.Halo(nout=187, base=wdir, load=True, is_gal=True)
hh = hmo.Halo(nout=187, base=wdir, load=True, is_gal=False)

good_gals = mma.close_gals(hh, gg, rscale=3)

# cluster 29176 ok_gals=[3, 5, 6, 9, 10, 11, 15, 16, 18, 20, 21, 22, 25, 26, 31, 32, 34, 37, 39, 43, 45, 49, 54, 55, 62, 63, 67, 68, 71]
#ok_gals =[2, 3,4, 6, 7, 10, 13, 18, 19, 22, 23, 28, 34,
#          40, 41, 42, 47, 48, 50, 51, 52, 54, 55, 59,
#          65, 66, 68, 73, 75, 81, 82, 83, 88, 89, 101, 102,
#          108, 113, 128, 131, 132, 135, 144, 148, 151, 157, 158,
#          170, 171, 178, 179,182,208, 214, 225, 421] # 10002
# Galaxies with good main progenitor list.
# But some galaxies have up-and-down mass evolution histories. 
# By looking at merger tree plots, I can roughly tell 
# major merger galaixies, minor merger galaxies, smooth accretion galaxies.


# In[23]:

print(" {:.2e}".format(min(alltrees.data['m'][alltrees.data['m'] > 0])))


# Load and compile catalogs

# In[6]:

alltrees = ctu.load_tree(wdir, is_gal=True)
# Serialize catalogs. -> Only main galaxies
# main galaxy list

ad = alltrees.data
tn = ad[ad['nout'] == nout_fi]

ad['r'] *= 2e2 # code unit to Mpc/h

cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout_fi) + '.pickle', 'rb'))
# Only good galaxies are saved in the catalogue.
# But there is no mass cut in good_gals list.
# And match_list_ind doest not handle cases where the array a is shorter than array b.
# Need another matching funciton that automatically excludes missing values.
good_gals = np.intersect1d(cat['id'], good_gals)
ind_ok = mtc.match_list_ind(cat['id'], good_gals, allow_swap=False)
cat = cat[ind_ok]

# Some galaxies are not in the tree. So can't follow their main progenitors.
cat = cat[(~np.isnan(cat['idx'])) * (cat['idx'] > 0)] 
# Some galaxies has no idx (-1). These are large gaalxies wihtout trees.
# I am keeping them hoping at some point they are useful in fixing bad trees.
idx_all = cat['idx'].astype(int) # fix so that idx is int from the beginning. 
#print(good_gals)


# In[7]:

mpgs = mma.compile_mpgs(alltrees, idx_all, wdir=wdir, cdir=cdir, nout_ini=nout_ini, nout_fi=nout_fi)


# Merger epochs

# In[9]:

mma.find_merger_epochs(alltrees, idx_all, mpgs, nout_ini=nout_ini, dist_gal_scale=1,
                      mass_ratio='early') 
# mass_ratio = early : satellite mass at when radii of two galaxies touch each other.
#            = max   : maximum mass of the satellite.   

for gal in mpgs:
    #print(gal.ids[0])
    filter_small_mergers(gal.merger)


# In[10]:

# Measure delta lambda
mma.measure_delta_lambda(mpgs, dt_before=7, dt_after=10, nout_ini=nout_ini)
mma.Maj_min_acc_ratio(mpgs)


# In[19]:

multipage_plots(mpgs, nout_ini=nout_ini, wdir=wdir, suffix=suffix)

pickle.dump(mpgs, open(dir_out + "mpgs" + suffix + ".pickle", 'wb'))

#plot_simple(mpgs, wdir=dir_out, suffix="")

plot_violin(mpgs, mstar_limit = 1e10)

plot_hist(mpgs, wdir=dir_out, suffix=suffix)


# ### Draw a vilon plot of Major/Minor/accretion contributio ratio

# ### Visual inspection. 
# Do I differentiate the effects of major/minor merger and accretion?
