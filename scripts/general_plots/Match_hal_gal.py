
# coding: utf-8

# # Moster plot for as many galaxies as possible

from tree import ctutils as ctu
from tree import halomodule as hmo
from load.info import Info
import pandas as pd
import pickle
import numpy as np
import matplotlib.pyplot as plt

# Load galaxy catalog
# and halo catalog
cluster = "29172"

wdir = "./" + cluster + '/' 
info = Info(187, base=wdir)


# In[6]:

def distv3d(halo, center):
    norm = np.sqrt(np.square(center['vx'] - halo['vx']) + 
                   np.square(center['vy'] - halo['vy']) + 
                   np.square(center['vz'] - halo['vz']))
    return norm

def distv(halo, center):
    norm = center['vx'] - halo['vx'] +            center['vy'] - halo['vy'] +            center['vz'] - halo['vz']
    return norm


def dist(halo, center):
    norm = np.sqrt(np.square(center['x'] - halo['x']) + 
                   np.square(center['y'] - halo['y']) + 
                   np.square(center['z'] - halo['z']))
    return norm 

def match_gal_hal_tree(gt, ht):
    nout = 187
    dt = 3 # compare progenitor at dt ago.
    
    gal_now = gt[gt["nout"]==nout]
    hal_now = ht[ht["nout"]==nout]
    
    gal_before = gt[gt["nout"]==nout-dt]
    hal_before = ht[ht["nout"]==nout-dt]    
    
    dominant = 0.1 # matched one should have less error by this amount or smaller 
                    # compared to the second best matched one.
    
    abs_tol_pos = 5e-5 # Position absolute tolerance [in code unit?]
    abs_tol_vel = 10   # velocity absolute tolerance [in kms?]
    
    for gal in gal_now:
        dd = dist(hal_now, gal)
        vv = distv(hal_now, gal)
        d_sort = np.argsort(dd)
        v_sort = np.argsort(vv)
        if (dd[d_sort[0]] < dominant * dd[d_sort[1]]) and (dd[d_sort[0]] < abs_tol_pos) and         (vv[v_sort[0]] < dominant * vv[v_sort[1]]) and (vv[v_sort[0]] < abs_tol_vel):
            gal['hosthalo'] = allhal.data['id'][d_sort[0]]
            i0.append(i)
            newhals[i] = allhal.data[d_sort[0]]
        else:
            atree = tree.atree(gt)
            prg = atree[dt]
            for gal2 in gal_before:
                dd = dist(hal_now, gal2)
                vv = distv(hal_now, gal2)
                d_sort = np.argsort(dd)
                v_sort = np.argsort(vv)
            
                

def get_comp_dist(hal_now, gal, nreturn=5):
    dd = dist(hal_now, gal)
    vv = distv(hal_now, gal)
    dd_q1 = np.percentile(dd,10)
    vv_q1 = np.percentile(vv,10)
    comp_dist = np.sqrt(np.square(dd/dd_q1) + np.square(vv/vv_q1))
    ind_sort = np.argsort(comp_dist)
    return comp_dist[ind_sort[:nreturn]], hal_now[ind_sort[:nreturn]]

def before_to_now(htdata, hals, dt):
    out = []
    for hal in hals:
        atree_hal = ctu.extract_main_tree_full(htdata, idx=hal['id'])
        out.append(atree_hal[hal['nout'] + dt])
    return out

def now_to_before(htdata, hals, dt):
    """
    progenitor of current halos.
    If does not exist, give -1
    """
    out =[]
    for hal in hals:
        idx = hal['id']
        try:
            #smalldata = htdata[htdata['tree_root_id'] == idx]
            #print(smalldata)
            atree_hal = ctu.extract_main_tree(htdata, idx=idx)
            out.append(atree_hal[dt])
        except:
            print("broken tree")
            out.append(-1)        
            
    return np.array(out)


# In[8]:

import tree.treemodule as tmo

from load.rd_GM import Gal
nout_fi = 187
#gt = pickle.load(open(wdir + "GalaxyMaker/Trees/extended_tree.pickle", "rb"))
#ht = pickle.load(open(wdir + "halo/Trees/extended_tree.pickle", "rb"))
gt = tmo.load_tree(wdir, is_gal=True)
ht = tmo.load_tree(wdir, is_gal=False)

dt = 3 # compare progenitor at dt ago.

gal_now = gt.data[gt.data["nout"]==nout_fi]
hal_now = ht.data[ht.data["nout"]==nout_fi]

#gal_before = gt.data[gt.data["nout"]==nout-dt]
hal_before = ht.data[ht.data["nout"]==nout_fi - dt]    

dominant = 0.4 # matched one should have less error by this amount or smaller 
                # compared to the second best matched one.

abs_tol_pos = 1e-2 # Position absolute tolerance [in Mpc]
abs_tol_vel = 100   # velocity absolute tolerance [in kms?]

comp_dists=[]
gal_ok = ctu.check_tree_complete(gt.data, nout_fi - dt, nout_fi, gal_now["id"], idx=True)

#comp_dists.append(comp_dist)

dts = [3,5]

hal_now = ht.data[ht.data["nout"]==nout_fi]

hal_3 = ht.data[ht.data["nout"]==nout_fi - 3]
hal_5 = ht.data[ht.data["nout"]==nout_fi - 5]
hal_this_list = [hal_3, hal_5]

result = []



i_gal_ok = []
for igal, gal in enumerate(gal_now):
    if gal['id'] not in gal_ok:
        continue
    else:
        i_gal_ok.append(igal)
        nreturn = 10
        comp_dist, good_hals_now = get_comp_dist(hal_now, gal, nreturn=10)
        # halo must be more massive than the galaxy
        matches=[]
        good_hals_now = good_hals_now[good_hals_now["m"] > gal["m"]]

        atree = ctu.extract_main_tree(gt.data, idx=gal['id'])

        matches.append(good_hals_now["Orig_halo_id"])
        for idt, dt in enumerate([3,5]):
            hal_this = hal_this_list[idt]
            gal_this = atree[dt]
            comp_dist_this, good_hals_this = get_comp_dist(hal_this, gal_this, nreturn=10)
            good_hals_this = good_hals_this[good_hals_this["m"] > gal_this["m"]]
            
            good_hals_prgsthis = now_to_before(ht.data, good_hals_now, dt)

            i_good = []
            i_good_prg=[]
            for i,ghthis in enumerate(good_hals_this['Orig_halo_id']):
                if ghthis in good_hals_prgsthis["Orig_halo_id"]:
                    #i_good_prg.append(i)
                    i_good.append(np.where(good_hals_prgsthis["Orig_halo_id"] == ghthis)[0][0])

            matches.append(good_hals_now["Orig_halo_id"][i_good])
        try:
            if matches[0][0] == matches[1][0] == matches[2][0]:
                matched = matches[0][0]
                #print(matched)
                result.append(matched)
            else:
                #print("Not good")
                result.append(-1)
        except:
            #print("Not good")
            result.append(-1)

        
result = np.array(result)


print( "Out of {} galaxies, failed to match {} galaxies.".format(len(result), sum(result < 0)))


#print(len(np.unique(result[result > 0])), sum(result > 0))
# There are multiple counted galaxies!


# load each galaxy
i_gal_ok = np.array(i_gal_ok)
good_gal = gal_now[i_gal_ok[result > 0]]["Orig_halo_id"]
good_hal = np.unique(result[result > 0])
#good_gal = [251,234]

for idgal in good_gal:
    gal = Gal(nout_fi, idgal, wdir=wdir, load=False)
    gal.load(cell="none", dm="none")
    # There are only two options: gm or raw. 
    # Others are ignored.
    
    gal.star['x'] -= gal.header['xg'][0]
    gal.star['y'] -= gal.header['xg'][1]
    gal.star['z'] -= gal.header['xg'][2]

    #xc,yc,zc = gm2code(gal.header['xg'], info)
    #gal.dm['x'] -= xc
    #gal.dm['y'] -= yc
    #gal.dm['z'] -= zc

    #gal.cell['x'] -= xc
    #gal.cell['y'] -= yc
    #gal.cell['z'] -= zc

    # rescale
    gal.star['x'] *= 1e3
    gal.star['y'] *= 1e3
    gal.star['z'] *= 1e3

    #gal.dm['x'] *= info.boxtokpc
    #gal.dm['y'] *= info.boxtokpc
    ##gal.dm['z'] *= info.boxtokpc
    
# measure stellar mass by  

# correct halo mass for extract substructure mass.


# In[ ]:

#gal.

