import matplotlib.pyplot as plt
import numpy as np
import load
import tree.halomodule as hmo
import pickle
import pandas as pd
import utils.match as mtc
import draw
from load.info import Info

def get_M1(z, M10=11.59 , M11=1.195):
    10**(M10 + M11(z/(z+1))


def m_M(M,z):
   N = get_M1
   N10 = 0.0351
   N11 = -0.0247

   beta = get_M1
   beta10 = 1.376
   beta11 = -0.826
 
   gamma = get_M1
   gamma10 = 0.608
   gamma11 = 0.329

   M1 = get_Ma(z, 11.59, 1.195)
   nn = N(z, N10, N11)
   bb = beta(z, beta10, beta11)
   gg = gamma(z, gamma10, gamma11)

   return 2 * N / ( (M/M1)**(-bb) + (M/M1)**(-gg))




xp = np.logspace(1e10, 1e15)
plt.plot(xp, moster_ref(xp))

clusters = ["05420", "39990", "01605", "05427", "36415",\
            "36413", "29176", "29172", "04466", "10002",\
            "17891", "24954", "28930", "32675", "35663",\
            "14172", "06098", "07206"]

fig, ax = plt.subplots(1,2)
fig.set_size_inches(18, 8)

mvir_all_c=[]
mstar_all_c=[]

info = Info(187, base='./05420/')
for cluster in clusters:
    wdir = "./" + cluster + '/'
    try:
        cat_final, ad_final = pickle.load(open(wdir + cluster + "moster.pickle", "rb"))
    except:
        continue

    # Color code central / satellite
    ind_cen = ad_final["level"] == 1
    ind_sat = ~ind_cen#np.where(ad_final["level"] != 1)[0]
    mstar = cat_final["mstar"]
    mvir = ad_final['mvir']

#    plt.clf()
    mvir_all_c.extend(mvir[ind_cen])
    mstar_all_c.extend(mstar[ind_cen])
    satellites = ax[0].scatter(np.log10(mvir[ind_sat]),\
                 mstar[ind_sat]/mvir[ind_sat] / (info.ob/info.om), \
                facecolors="blue", edgecolors="blue",\
                label="satellite", )
    centrals = ax[0].scatter(np.log10(mvir[ind_cen]),\
                 mstar[ind_cen]/mvir[ind_cen] / (info.ob/info.om),\
                facecolors = "red", edgecolors="red",\
                label="central")
    
    centrals = ax[1].scatter(np.log10(mvir[ind_cen]),\
                 mstar[ind_cen]/mvir[ind_cen] / (info.ob/info.om),\
                facecolors = "red", edgecolors="red",\
                label="central")

ax[0].set_ylim([0,5])
ax[1].set_ylim([0,2])

#    ax[1].set_title(cluster)
#    ax = plt.gca()

#plt.legend(handles=[centrals, satellites])
#plt.title(cluster)
ax[0].set_ylabel(r"$ M_{\star} / M_{200} / (\Omega_{b} / \Omega_{m} )$")
ax[0].set_xlabel(r"log$[M_{200} / M_{\odot}]$")

ax[1].set_ylabel(r"$ M_{\star} / M_{200} / (\Omega_{b} / \Omega_{m} )$")
ax[1].set_xlabel(r"log$[M_{200} / M_{\odot}]$")

plt.savefig("./ALL_Moster_plot.png")
#plt.savefig(wdir + cluster + "Moster_plot_central_only.png")

plt.close()
#fig, ax = plt.subplots(1,2)
#plt.clf()
plt.scatter(np.log10(mvir_all_c), np.log10(mstar_all_c))
plt.savefig("./Mstar_Mhal_cen.png")
