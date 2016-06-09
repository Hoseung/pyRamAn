# -*- coding: utf-8 -*-
"""
Created on Wed May 18 23:52:00 2016

@author: hoseung
"""

import numpy as np
import pickle
import analysis.evol_lambda as evl
import tree.ctutils as ctu
import pandas as pd


def load_cat(fname):
    return pd.DataFrame(pickle.load(open(fname, "rb"))).to_records()


class Merger():
    pass


def get_tick_locations(org_range, x_range):
    nticks = 5
    xx = x_range.ptp()
    dt = xx // (nticks - 1)
    x_ticks = np.ceil(x_range[0] + dt * np.arange(nticks))
    x_ticks = x_ticks[x_ticks < max(x_range)]
    labels = x_ticks[x_ticks < max(x_range)]
    
    tick_pos = (labels - x_range[0]) * org_range.ptp()/xx + org_range[0]
    return  tick_pos, labels



def find_merger_epochs(alltrees, idx_all, mpgs, nout_ini=37, dist_gal_scale=2,
                       mass_ratio="early"):
    """
    
        parameter
        ---------
        dist_gal_scale 
            if two galaxies are closer than dist_gal_scale * (sum of raidus of the two),
            that epoch is the nout_init_merger.
        nout_ini
            blabla
    """
    verbose=False
    gal_list=[]
    mr_list=[]
    nout_list=[]
    nout_ini_list=[] # initial time when two halos(Galaxy stellar components in this case) overlap. 

    for idx in idx_all:
        # full tree of a galaxy
        atree = ctu.extract_a_tree(alltrees.data, idx)

        # main progenitor tree
        main = ctu.extract_main_tree(atree, idx)

        x_nout = main['nout'].flatten()
        i_nout_ok = x_nout > nout_ini
        main = main[i_nout_ok]
        #x_nout = x_nout[i_nout_ok]
        pos = np.zeros((3,len(main)))
        pos[0,:] = main['x']
        pos[1,:] = main['y']
        pos[2,:] = main['z']


        mass_ratios_single = np.zeros(len(main))
        nout_inits = np.zeros(len(main))
        for i, nout in enumerate(main['nout']):
            # merger ratio
            i_prgs = np.where(atree['desc_id'] == main['id'][i])[0]

            # multiple prgs = merger
            if len(i_prgs) > 1:
                if verbose: print(" {} mergers at nout = {}".format(len(i_prgs), nout))
                id_prgs = atree['id'][i_prgs]
                mass_prgs = atree['m'][i_prgs]
                m_r = mass_prgs / max(mass_prgs)

                i_sat = np.argmax(mass_prgs[1:]) + 1
                #mass_ratios_single[i] = max([mass_prgs[1:] / max(mass_prgs)][0])
                if mass_ratio == "max":
                    mass_ratios_single[i] = mass_prgs[i_sat] / mass_prgs[0]

                satellite = ctu.extract_main_tree(atree, id_prgs[i_sat], no_subset=True)

                nout_min = max([min(main['nout']), min(satellite['nout'])])
                i_main_ok = (main['nout'] > nout_min) * (main['nout'] < nout)
                i_sat_ok = (satellite['nout'] > nout_min) * (satellite['nout'] < nout)
                satellite = satellite[i_sat_ok]

                dd = np.sqrt(np.square(pos[0,i_main_ok] - satellite['x']) \
                             + np.square(pos[1,i_main_ok] - satellite['y'])
                             + np.square(pos[2,i_main_ok] - satellite['z']))
                rgal_tot = main['r'][i_main_ok] + satellite['r']
                #print(" Galaxy sizes : main {}, and the second {}, and the sum {}".format(
                #        main['r'][i_main_ok], satellite['r'], rgal_tot))
                #print(" dd :", dd)
                if sum(dist_gal_scale * rgal_tot < dd) > 0:
                    #print(50 * rgal_tot - dd)
                    #print(satellite['nout'][50 * rgal_tot < dd])
                    nout_inits[i] = max(satellite['nout'][dist_gal_scale * rgal_tot < dd])
                    if mass_ratio == "early":
                        
                        mass_ratios_single[i] = satellite['m'][satellite['nout'] == nout_inits[i]] / mass_prgs[0]
                        
                else:
                    nout_inits[i] = nout
                if verbose:
                    print(" Mass ratios : ", m_r, nout, nout_inits[i])

            else:
                mass_ratios_single[i] = 0

        ind_ok = np.where(mass_ratios_single > 0.01)[0]
        if len(ind_ok) > 0:
            # if a satellite oscillates around the host, 
            # it could be identified as multiple mergers with short time interval. 
            # leave only the first passage / merger.
            # No, it doesn't happen in ConsistentTrees.

            #good =[]
            #for i in range(len(ind_ok)-1):
            #    if ind_ok[i+1] > ind_ok[i] + 2:
            #        good.append(ind_ok[i])
            #good.append(ind_ok[-1])
            #ind_ok = good
            mr = 1./mass_ratios_single[ind_ok]

            gal_list.append(idx)
            mr_list.append(mr)
            nout_list.append(x_nout[ind_ok])    
            nout_ini_list.append(nout_inits[ind_ok])


    inds=[]
    for i, gal in enumerate(mpgs):
        galid = gal.data['idx'][0]
        ind = np.where(galid == gal_list)[0]
        if len(ind) > 0:
            inds.append(i)
            merger = Merger()
            merger.mr = mr_list[ind]
            merger.nout = nout_list[ind]
            merger.nout_ini = nout_ini_list[ind]
            gal.merger = merger
        else:
            gal.merger = None
    #return gal_list, mr_list, nout_list, nout_ini_list

    

def close_gals(halo, gals, return_ind=True, rscale=3.0):
    import numpy as np
    i_cluster = np.argmax(halo.data['np'])
    cluster = halo.data[i_cluster]
    xc = cluster['x']
    yc = cluster['y']
    zc = cluster['z']

    dd = np.square(gals.data['x'] - xc) + \
         np.square(gals.data['y'] - yc) + \
         np.square(gals.data['z'] - zc)
    
    r_clu = rscale * cluster['rvir']
    #print(cluster['rvir'])
    if return_ind:
        return np.where(dd < r_clu**2)[0]
    else:
        return gals.data['id'][dd < r_clu**2][0]
    
    

def compile_mpgs(alltrees, idx_all, wdir='./', cdir='easy/', nout_ini=37, nout_fi=187):
    ad = alltrees.data
    mpg_tmp = []
    for i, idx in enumerate(idx_all):
        mpg_tmp.append(evl.MainPrg(ad, idx))
        #print(i, idx)

    for nout in range(nout_ini, nout_fi + 1):
        cat = load_cat(wdir + cdir + 'catalog' + str(nout) + '.pickle')
#        cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout) + '.pickle', 'rb'))
        for gal in mpg_tmp:
            gal.set_data(cat, nout)
        #print(nout)

    mpgs = []
    while len(mpg_tmp) > 0:
        mpgs.append(mpg_tmp.pop())
        
    return mpgs



def Maj_min_acc_ratio(mpgs, dt=5, major_ratio=3):
    for igal, gal in enumerate(mpgs):
        delta_lambda_tot = np.average(gal.data['lambda_r'][:dt]) - np.average(gal.data['lambda_r'][-dt:])

        delta_lambda_major=0
        delta_lambda_minor=0
        if gal.merger is not None:

    #        for mr, xx, delta in zip(gal.merger.mr, gal.merger.nout, gal.merger.delta):
    #            ax[0].scatter(mr, delta)
            #print(gal.merger.mr, gal.merger.delta)
            i_major = np.where(gal.merger.mr <= major_ratio)[0]
            if len(i_major) > 0:
                #print(gal.merger.delta, i_major, "M")
                #print(gal.merger.delta[0])
                delta_lambda_major = np.sum(gal.merger.delta[i_major])

            i_minor = gal.merger.mr > major_ratio
            if sum(i_minor) > 0:
                delta_lambda_minor = sum(gal.merger.delta[i_minor])
                #print(gal.merger.delta, i_minor, "m")


        delta_lambda_other = delta_lambda_tot - delta_lambda_major - delta_lambda_minor
        gal.dlt = delta_lambda_tot
        gal.dlo = delta_lambda_other
        gal.dlM = delta_lambda_major
        gal.dlm = delta_lambda_minor

  
        
        
def measure_delta_lambda(mpgs, dt_before=7, dt_after=7, nout_ini=37,
                         filter_smaller=True):
    """
        Note that nout are in descending order.
        
        physical meaning of dt_after is the time for a merger to settle down.
        And, if there are multiple mergers within dt_after 
        the effect of the larger merger and the rest are all mixed up. 
        Then I assume that it is more reasonable to take the 'mixed' effect
        is of the larger merger alone.
        
    """
    from scipy.signal import medfilt
    for gal in mpgs:
        ind_nout = gal.nouts > nout_ini
        gal.nouts = gal.nouts[ind_nout]
        gal.data = gal.data[ind_nout]
        gal.smoothed_lambda = medfilt(gal.data['lambda_r'], kernel_size=5)

        if gal.merger is not None:
            #if filter_smaller:            
                #gal.merger = filter_smaller_mergers(gal.merger)
            delta_lambda =[]
            
            for mr, xx, x2 in zip(gal.merger.mr, gal.merger.nout, gal.merger.nout_ini):
                i_nout = np.where(gal.nouts == xx)[0]
                iini_nout = np.where(gal.nouts == x2)[0]

                lambda_after = np.average(gal.smoothed_lambda[max([0, i_nout - dt_after]) : i_nout])
                lambda_before = np.average(gal.smoothed_lambda[iini_nout:min([len(gal.data), iini_nout + dt_before])])
                delta_lambda.append(lambda_after - lambda_before)
            gal.merger.delta = np.array(delta_lambda)
