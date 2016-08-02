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

class Merger():
    pass


def load_cat(fname):
    return pd.DataFrame(pickle.load(open(fname, "rb"))).to_records()


def get_tick_locations(org_range, x_range):
    nticks = 5
    xx = x_range.ptp()
    dt = xx // (nticks - 1)
    x_ticks = np.ceil(x_range[0] + dt * np.arange(nticks))
    x_ticks = x_ticks[x_ticks < max(x_range)]
    labels = x_ticks[x_ticks < max(x_range)]
    
    tick_pos = (labels - x_range[0]) * org_range.ptp()/xx + org_range[0]
    return  tick_pos, labels



def find_merger_epochs(alltrees, 
                        idx_all, 
                        mpgs, 
                        nout_ini=37, 
                        dist_gal_scale=2,
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
                             + np.square(pos[2,i_main_ok] - satellite['z'])) * 1e3 # in kpc
                rgal_tot = main['rvir'][i_main_ok] + satellite['rvir']
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
            #print(gal.merger.mr, gal.merger.delta_l)
            i_major = np.where(gal.merger.mr <= major_ratio)[0]
            if len(i_major) > 0:
                #print(gal.merger.delta, i_major, "M")
                #print(gal.merger.delta[0])
                delta_lambda_major = np.sum(gal.merger.delta_l[i_major])

            i_minor = gal.merger.mr > major_ratio
            if sum(i_minor) > 0:
                delta_lambda_minor = sum(gal.merger.delta_l[i_minor])
                #print(gal.merger.delta, i_minor, "m")


        delta_lambda_other = delta_lambda_tot - delta_lambda_major - delta_lambda_minor
        gal.dlt = delta_lambda_tot
        gal.dlo = delta_lambda_other
        gal.dlM = delta_lambda_major
        gal.dlm = delta_lambda_minor



def multipage_plot(mpgs, nout_ini=37, nout_fi = 187,
                    wdir='./',
                    suffix="",
                    dt_after = 1.0,
                    dt_before = 1.0,
                    do_plot=True):

    
    from matplotlib.backends.backend_pdf import PdfPages
    fig, ax = plt.subplots(1,2, sharex=True)
    plt.subplots_adjust(hspace=0.001)
    fig.set_size_inches(8,4)

    with PdfPages(wdir + 'MMA_plots' + suffix +'.pdf') as pdf:
        measure_delta(mpgs, nout_ini=37, nout_fi = 187,
                    wdir='./',
                    dt_after = dt_after,
                    dt_before = dt_before,
                    ax=ax)
        
        
        pdf.savefig()
        ax[0].clear()
        ax[1].clear()

    #plt.show()
    plt.close()
    



def measure_delta(mpgs, nout_ini=37, nout_fi = 187,
                    wdir='./',
                    dt_after = 1.0,
                    dt_before = 1.0,
                    ax=None):
    """
    Measure lambda change at every merger and draw lambda evolution with merger events.
    time span over which average lambda value is measured is given in Gyr unit.


    smoothing is done elsewhere.         
    """

    from scipy.signal import medfilt
    from utils.util import dgyr2dnout    
    for gal in mpgs:
        ind_nout = gal.nouts > nout_ini
        gal.nouts = gal.nouts[ind_nout]
        gal.data = gal.data[ind_nout]
#        gal.smoothed_lambda = medfilt(gal.data['lambda_r'], kernel_size=5)
        if ax is not None:
            ax[0].plot(gal.nouts, np.log10(gal.data['mstar']), 'b:')
            ax[0].set_xlim([50,190])
            ax[0].set_title(str(gal.ids[0]) + ", " + str(gal.idxs[0]))
            ax[1].plot(gal.nouts, gal.smoothed_lambda, 'black')

        if gal.merger is not None:
            delta_lambda =[]
            delta_mass = []

            for mr, xx, x2 in zip(gal.merger.mr, gal.merger.nout, gal.merger.nout_ini):
                # index of merger (final coalescence)
                i_nout = np.where(gal.nouts == xx)[0]
                # index of the begining of merger (as the time two galaxies overlap
                # => rga1 + rgal2 < distance)
                iini_nout = np.where(gal.nouts == x2)[0]
                # snapshot number away by dt from nout
                nout_after = dgyr2dnout(dt_after, xx)
                nout_before = dgyr2dnout(-1 * dt_before, x2)
                # indices of the snapshots.
                inout_after = np.where(gal.nouts == nout_after)[0]
                inout_before = np.where(gal.nouts == nout_before)[0]

                # if time_{merger} + dt > t_{z=0}
                if i_nout == 0:
                    l_inds_after = [0]
                else:
                    l_inds_after = range(max([0, inout_after]), i_nout)

                nouts_after = gal.nouts[l_inds_after]
                l_after = gal.smoothed_lambda[l_inds_after]
                m_after = gal.data['mstar'][l_inds_after]
                lambda_after = np.average(l_after)
                mass_after = np.average(m_after)


                l_inds_before = range(iini_nout,min([len(gal.data), inout_before]))
                nouts_before = gal.nouts[l_inds_before]
                l_before = gal.smoothed_lambda[l_inds_before]
                m_before = gal.data['mstar'][l_inds_before]
                lambda_before = np.average(l_before)
                mass_before = np.average(m_before)

                delta_lambda.append(lambda_after - lambda_before)
                delta_mass.append(mass_after - mass_before)


                if ax is not None:
                    ax[1].plot(nouts_after,l_after, 'g-')                                           
                    # Check again.
                    nn = range(min(nouts_after) - 5, min([nout_fi, max(nouts_after) + 5]))
                    ax[1].plot(nn, [lambda_after]*len(nn), "g:")                    

                    ax[1].plot(nouts_before,l_before, 'r-')
                    nn = range(max([nout_ini, min(nouts_before) - 5]), max(nouts_before) + 5)
                    ax[1].plot(nn, [lambda_before]*len(nn), "r:")                    

                    ax[0].axvline(xx, linestyle=':')
                    ax[0].annotate("{:.1f}".format(mr), xy=(xx,0.8))
                    ax[1].axvline(xx, linestyle=':')
                    ax[1].annotate("{:.1f}".format(mr), xy=(xx,0.8))
                    ax[1].axvline(xx, linestyle=':')
                    ax[1].axvline(x2, linestyle=':', c='g')

            gal.merger.delta_l = np.array(delta_lambda)
            gal.merger.delta_m = np.array(delta_mass)        



def smooth(x, beta=5, window_len=20, monotonic=False, clip_tail_zeros=True):
    """ 
    kaiser window smoothing.
    
    If len(x) < window_len, window_len is overwritten to be len(x).
    This ensures to return valid length fo array, but with modified window size.
       
    
    beta = 5 : Similar to a Hamming
    
    
    """
    if clip_tail_zeros:
        x = x[:max(np.where(x > 0)[0])+1]
    
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
    window_len = min([window_len, len(x)])
    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]] # concatenate along 0-th axis.
    # periodic boundary.
    w = np.kaiser(window_len,beta)
    y = np.convolve(w/w.sum(), s, mode='valid')
    if monotonic: 
        return y[int(window_len)/2:len(y)-int(window_len/2) + 1] + xx
    else:
        return y[int(window_len)/2:len(y)-int(window_len/2) + 1]
