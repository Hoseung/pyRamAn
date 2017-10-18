
# coding: utf-8

# # Moster plot for as many galaxies as possible

# In[135]:

class Dummy():
    def __init__(self):
        pass

def radial_profile_cut(star,
                       den_lim=2e6, den_lim2=5e6,
                       mag_lim=25, nbins=100, rmax=50, dr=0.5,
                       debug=False):
    # 2D photometry. (if rotated towards +y, then use x and z)
    # now assuming +z alignment. 
    xx = star['x']
    yy = star['y']
    zz = star['z']
    vx = star['vx']
    vy = star['vy']
    vz = star['vz']
    mm = star['m']
    
    meta = Dummy()
    
    rr = np.sqrt(np.square(xx) + np.square(yy))# in kpc unit
    if debug:
        print(min(rr), max(rr), min(xx), max(xx))

    # Account for weights.
    i_sort = np.argsort(rr)
    r_sorted = rr[i_sort]
    m_sorted = mm[i_sort]

    rmax = min([np.max(rr), rmax])
    nbins = int(rmax/dr)

    frequency, bins = np.histogram(r_sorted, bins = nbins, range=[0, rmax])
    bin_centers = bins[:-1] + 0.5 * dr # remove the rightmost boundary.

    m_radial = np.zeros(nbins)
    ibins = np.concatenate((np.zeros(1), np.cumsum(frequency)))

    i_r_cut1 = nbins -1 # Maximum value
    # on rare occasions, a galaxy's stellar surface density
    # never crosses the density limit. Then i_r_cut1 = last index.
    for i in range(nbins):
        m_radial[i] = np.sum(m_sorted[ibins[i]:ibins[i+1]])
        if (m_radial[i]/(2 * np.pi * bin_centers[i] * dr)) < den_lim:
            i_r_cut1 = i-1
            break
    #i_r_cut2= np.argmax(m_radial/(2 * np.pi * bin_centers * dr) < den_lim2)
    # If for some reason central region is less dense,
    # profile can end at the first index.
    # Instead coming from backward, search for the point the opposite condition satisfied.
    den_radial_inverse = m_radial[::-1]/(2 * np.pi * bin_centers[::-1] * dr)
    
    if max(den_radial_inverse) < 1.5 * den_lim2:
        #print("Not dense enough")
        return False
    i_r_cut2=len(m_radial) - np.argmax(den_radial_inverse > den_lim2) -1
    if debug:
        print("[galaxy.Galaxy.radial_profile_cut] m_radial \n", m_radial)
        print("[galaxy.Galaxy.radial_profile_cut] den_radial_inverse \n", den_radial_inverse)
        print("[galaxy.Galaxy.radial_profile_cut] i_r_cut2", i_r_cut2)

    mtot2 = sum(m_radial[:i_r_cut2])
    mtot1 = sum(m_radial[:i_r_cut1])
    i_reff2 = np.argmax(np.cumsum(m_sorted) > (0.5*mtot2))
    i_reff1 = np.argmax(np.cumsum(m_sorted) > (0.5*mtot1))
    meta.reff2 = r_sorted[i_reff2]
    meta.reff  = r_sorted[i_reff1]
    #print(bin_centers, i_r_cut2, m_radial)
    meta.rgal2 = max([bin_centers[i_r_cut2],4*meta.reff2])
    meta.rgal  = max([bin_centers[i_r_cut1],4*meta.reff])#bin_centers[i_r_cut1]

    if debug: print("[galaxy.Galaxy.radial_profile_cut] mtot, mtot2", mtot1, mtot2)

    i_close = i_sort[:np.argmax(np.cumsum(m_sorted) > (0.2*mtot2))] # 20% closest particles
#        i_close = np.argsort(rr)[:min([i_reff1])]
#        i_close = i_sort[:min([i_reff1])]

    meta.mtot = mtot2

    meta.vxc = np.average(vx[i_close])
    meta.vyc = np.average(vy[i_close])
    meta.vzc = np.average(vz[i_close])

    return meta



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


# In[133]:

def measure_mstar_mhal(cluster, 
                       dominant = 0.4,
                       abs_tol_pos = 1e-2,
                       abs_tol_vel = 100,
                       nout_fi = 187,
                       dts = [3,5]):
    wdir = "./" + cluster + '/' 
    info = Info(187, base=wdir)

    gt = tmo.load_tree(wdir, is_gal=True)
    ht = tmo.load_tree(wdir, is_gal=False)

    gal_now = gt.data[gt.data["nout"]==nout_fi]
    hal_now = ht.data[ht.data["nout"]==nout_fi]

    comp_dists=[]
    gal_ok = ctu.check_tree_complete(gt.data, nout_fi - max(dts), nout_fi, gal_now["id"], idx=True)

    #comp_dists.append(comp_dist)

    hal_now = ht.data[ht.data["nout"]==nout_fi]

    hal_3 = ht.data[ht.data["nout"]==nout_fi - 3]
    hal_5 = ht.data[ht.data["nout"]==nout_fi - 5]
    hal_this_list = [hal_3, hal_5]

    result = []
    mhal_result = []
    dist_error = []

    i_gal_ok = []
    ok_gals = []
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
                    #result.append(matched)
                    result.append(good_hals_now[0])
                    ok_gals.append(gal)
                    #mhal_result.append(good_hals_now["mvir"][0])
                    #dist_error.append(comp_dist_this[0])
                #else:
                    #print("Not good")
                    #result.append(-1)
                    #mhal_result.append(-1)
                    #dist_error.append(-1)
            except:
                pass
                #print("Not good")
                #result.append(-1)
                #mhal_result.append(-1)
                #dist_error.append(-1)


    result = np.array(result)
    ok_gals = np.array(ok_gals)
    #mhal_result = np.array(mhal_result)
    #dist_error = np.array(dist_error)

    print( "Out of {} galaxies, matched {} galaxies.".format(len(gal_now), len(result)))

    # filter duplicates
    unq, unq_idx, unq_cnt = np.unique(result["Orig_halo_id"], return_inverse=True, return_counts=True)
    cnt_mask = unq_cnt > 1
    cnt_idx, = np.nonzero(cnt_mask)
    idx_mask = np.in1d(unq_idx, cnt_idx)
    idx_idx, = np.nonzero(idx_mask)
    srt_idx = np.argsort(unq_idx[idx_mask])
    dup_idx = np.split(idx_idx[srt_idx], np.cumsum(unq_cnt[cnt_mask])[:-1])

    # Remove smaller duplicates and leave the largest galaxy.
    remove_inds=[]
    for dup in dup_idx:
        ind_all = np.full(len(dup), True, dtype=bool)
        #print(ind_all)
        imax = np.argmax(ok_gals["m"][dup])
        ind_all[imax] = False
        remove_inds.extend(dup[ind_all])

    remove_inds = np.array(remove_inds)
    inds_ok = np.full(len(result), True, dtype=bool)
    inds_ok[remove_inds] = False


    # load each galaxy and measure stellar mass
    good_gal = ok_gals["Orig_halo_id"][inds_ok]
    mhal_result = result["mvir"][inds_ok]

    mstar = []
    mhal = []

    for mhal_this, idgal in zip(mhal_result, good_gal):
        gal = Gal(nout_fi, idgal, wdir=wdir, load=False)
        gal.load(cell="none", dm="none")
        # There are only two options: gm or raw. 
        # Others are ignored.

        gal.star['x'] -= gal.header['xg'][0]
        gal.star['y'] -= gal.header['xg'][1]
        gal.star['z'] -= gal.header['xg'][2]

        # rescale
        gal.star['x'] *= 1e3
        gal.star['y'] *= 1e3
        gal.star['z'] *= 1e3
        gal.star['m'] *=1e11
        gal.meta = radial_profile_cut(gal.star)
        if gal.meta is not False:
            mstar.append(gal.meta.mtot)
            mhal.append(mhal_this)
        else:
            pass
    
    return mstar, mhal



    # correct halo mass for extract substructure mass.


# In[134]:

from tree import ctutils as ctu

from tree import halomodule as hmo
from load.info import Info
import pandas as pd
import pickle
import numpy as np
import matplotlib
matplotlib.use("Qt4Agg")
import matplotlib.pyplot as plt
import tree.treemodule as tmo
from load.rd_GM import Gal

# Load galaxy catalog
# and halo catalog
clusters = ["29172"]
#nout_fi = 187

#dominant = 0.4 # matched one should have less error by this amount or smaller 
                # compared to the second best matched one.

#abs_tol_pos = 1e-2 # Position absolute tolerance [in Mpc]
#abs_tol_vel = 100   # velocity absolute tolerance [in kms?]

#dts = [3,5]

mstar_all = []
mhal_all = []

for cluster in clusters:
    mstar_t, mhal_t = measure_mstar_mhal(cluster, 
                                           dominant = 0.4,
                                           abs_tol_pos = 1e-2,
                                           abs_tol_vel = 100,
                                           nout_fi = 187,
                                           dts = dts)
    mstar_all.extend(mstar_t)
    mhal_all.extend(mhal_t)

    
    


# In[ ]:



