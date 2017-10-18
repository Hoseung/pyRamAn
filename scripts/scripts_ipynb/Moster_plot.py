
# coding: utf-8

# # M* vs Mhalo

# Again, galaxy - halo matching is required.

# In[18]:

class Dummy():
    def __init__(self):
        pass


def region_from_xyz(x,y,z, r=0):
    """
    create a region to include all points.
    
    Parameters
    ----------
    x: float array
    y: float array
    z: float array
    r: float array, optional
        If given, the region encloses all points as spheres with radius r.
    """
    xmi = min(x-0)
    xma = max(x+0)
    ymi = min(y-0)
    yma = max(y+0)
    zmi = min(z-0)
    zma = max(z+0)
    import utils.sampling as smp
    return smp.set_region(ranges=[[xmi, xma],[ymi,yma], [zmi,zma]])


def get_tick_points(xmin, xmax):
    dx = xmax - xmin
    dxx = round(dx/5)

    mid = round((xmin + xmax)/2)
    xstart = mid - 2*dxx
 
    ticks = np.arange(xstart, xmax, dxx)
    pos = ((ticks - xmin)/dx )
   
    return pos, ticks


def ticks_from_region(ranges, ax, npix):
    """
    modify ticks of given ax according to region.
    xmin, xmax, ymin, ymax must be set according to the region.
    """
    xmin = ranges[0][0]
    xmax = ranges[0][1]
    ymin = ranges[1][0]
    ymax = ranges[1][1]
    
    xpos, xticks = get_tick_points(xmin, xmax)
    ypos, yticks = get_tick_points(ymin, ymax)
    xpos *= npix
    ypos *= npix
       
    ax.set_xticks(xpos)
    ax.set_yticks(ypos)
    ax.set_xticklabels([str(i) for i in xticks])
    ax.set_yticklabels([str(i) for i in yticks])


# In[22]:

def distv(halo, center):
    norm = np.sqrt(np.square(center['vx'] - halo.vx) + 
                   np.square(center['vy'] - halo.vy) + 
                   np.square(center['vz'] - halo.vz))
    return norm

def dist(halo, center):
    norm = np.sqrt(np.square(center['x'] - halo.x) + 
                   np.square(center['y'] - halo.y) + 
                   np.square(center['z'] - halo.z))
    return norm    


def associate_gal_hal(allgal, allhal):
    # associate halos with galaxies. 
    i0=[]
    i1=[]
    # equivalent with allhal.data
    newhals = np.recarray(len(allgal.data), dtype=allhal.data.dtype) 
    
    for i, gal in enumerate(allgal.data):
        dd = dist(allhal.data, gal)
        d_sort = np.argsort(dd)
        if (dd[d_sort[0]] < 0.1 * dd[d_sort[1]]) and (dd[d_sort[0]] < 5e-5):
            gal['hosthalo'] = allhal.data['id'][d_sort[0]]
            i0.append(i)
            newhals[i] = allhal.data[d_sort[0]]
        elif (dd[d_sort[0]] < 0.3 * dd[d_sort[1]]) and (dd[d_sort[0]] < 2.5e-5):
            gal['hosthalo'] = allhal.data['id'][d_sort[0]]
            i0.append(i)
            newhals[i] = allhal.data[d_sort[0]]
        else:
            dv = distv(allhal.data, gal)
            d_nomi = dd < 2e-4 # within 40kpc!! 
            v_nomi = dv < 150
            
            #print(i, len(dd), len(dv))
            #if len(d_nomi) != len(v_nomi):
            #    continue
            #print(d_nomi, v_nomi)
            #print(len(d_nomi), len(v_nomi))
            ok = d_nomi * v_nomi
            if sum(ok) > 0:
                dddv = dd[ok]/2e-4 + dv[ok]/150
                dddv_sort = np.argsort(dddv)
                gal['hosthalo'] = allhal.data['id'][dddv_sort[0]]
                i0.append(i)
                newhals[i] = allhal.data[dddv_sort[0]]
            else:
                gal['hosthalo'] = -1
                i1.append(i)
                newhals[i]['x'] = gal['x']
                newhals[i]['y'] = gal['y']
                newhals[i]['z'] = gal['z']
                newhals[i]['vx'] = gal['vx']
                newhals[i]['vy'] = gal['vy']
                newhals[i]['vz'] = gal['vz']
                newhals[i]['r'] = gal['r'] * 2 # arbitrary
                newhals[i]['m'] = gal['m'] * 2 # arbitrary
                newhals[i]['rvir'] = gal['rvir'] * 2 # arbitrary
                newhals[i]['mvir'] = gal['mvir'] * 2 # arbitrary
                # small radius assumes most of the DM, gas have been lost.
                # star is the only major component.
                newhals[i]['id'] = -1 * gal['id'] # negative ID = fake halo.
    
    allhal.data = newhals
    return allhal


# In[170]:

def plot_gal_hal(gal, hal, ind1=None, ind2=None, rscale=1.0, npix=200):
    """
    hal = gal associated halo catalog (matchted).
    """
    fig, ax = plt.subplots()
    fig.set_size_inches(16,16)
    md = hal.data 
    region = region_from_xyz(md['x'], md['y'], md['z'], md['r'])
    draw.pp.pp_halo(hal, ind=ind1, region=region, npix=npix, rscale=rscale, axes=ax, colors=500)
    ranges = np.array(region['ranges']) * info.pboxsize
    ticks_from_region(ranges, ax, npix)
    draw.pp.pp_halo(gal, ind=ind2, region=region, npix=npix, rscale=rscale, axes=ax, colors=10)
    plt.savefig(wdir + cluster + "halo_galaxy_match.png")
    ax.set_aspect('equal')
    plt.close()
    


# In[182]:

#import matplotlib
#matplotlib.use("Qt4Agg")
import matplotlib.pyplot as plt
import numpy as np
import load
import tree.halomodule as hmo
import pickle
import pandas as pd
import utils.match as mtc
import draw
from load.info import Info

clusters = ["05420", "39990", "01605", "05427", "36415",\
            "36413", "29176", "29172", "04466", "10002",\
            "17891", "24954", "28930", "32675", "35663",\
            "14172", "06098", "07206"] 
#cluster = clusters[0]
for cluster in clusters:
    wdir = "./" + cluster + '/'
    info = Info(187, base=wdir)
    print(cluster)
    try:
        cat = pd.DataFrame(pickle.load(open(wdir         + "./easy_new/catalog187.pickle", "rb"))).to_records()
    except:
        print("No lambda catalog", cluster)
        continue

    hh = hmo.Halo(187, base=wdir, is_gal=False)
    gg = hmo.Halo(187, base=wdir, is_gal=True)
    # match galaxy and halo
    allhal = associate_gal_hal(gg,hh)
    ad = allhal.data
    gd = gg.data
    # Match map
    plot_gal_hal(gg, allhal)


    # Galaxies in the lambda catalog
    ind_g = mtc.match_list_ind(gd["id"], cat["id"])
    final_gals = gd["id"][ind_g]
    # filter out failed ones.
    good = (cat['id'] == final_gals) * cat['mstar'] > 0

    # host halos.
    ind_d = mtc.match_list_ind(ad["id"], gd["hosthalo"][ind_g[good]])

    # Final samples
    cat_final = cat[good]
    ad_final = ad[ind_d]

    # Color code central / satellite
    ind_cen = np.where(ad_final["level"] == 1)[0]
    ind_sat = np.where(ad_final["level"] != 1)[0]
    mstar = cat_final["mstar"]
    mvir = ad_final['mvir']

    plt.clf()
    satellites = plt.scatter(np.log10(mvir[ind_sat]),\
                 mstar[ind_sat]/mvir[ind_sat] / (info.ob/info.om), \
                facecolors="blue", edgecolors="blue",\
                label="satellite")
    centrals = plt.scatter(np.log10(mvir[ind_cen]),\
                 mstar[ind_cen]/mvir[ind_cen] / (info.ob/info.om),\
                facecolors = "red", edgecolors="red",\
                label="central")
    plt.legend(handles=[centrals, satellites])
    plt.title(cluster)
    ax = plt.gca()
    ax.set_ylabel(r"$ M_{\star} / M_{200} / (\Omega_{b} / \Omega_{m} )$")
    ax.set_xlabel(r"log$[M_{200} / M_{\odot}]$")
    plt.savefig(wdir + cluster + "Moster_plot.png")
    
    plt.clf()
    centrals = plt.scatter(np.log10(mvir[ind_cen]),\
                 mstar[ind_cen]/mvir[ind_cen] / (info.ob/info.om),\
                facecolors = "red", edgecolors="red",\
                label="central")
    plt.title(cluster)
    ax = plt.gca()
    ax.set_ylabel(r"$ M_{\star} / M_{200} / (\Omega_{b} / \Omega_{m} )$")
    ax.set_xlabel(r"log$[M_{200} / M_{\odot}]$")
    plt.savefig(wdir + cluster + "Moster_plot_central_only.png")

    pickle.dump([cat_final, ad_final], open(wdir + cluster + "moster.pickle", "wb"))

# In[177]:




# In[179]:



# In[ ]:



