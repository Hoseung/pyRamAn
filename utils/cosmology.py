# -*- coding: utf-8 -*-
"""
cosmology utils.

... use astropy.cosmology. that is a full furnished util. 

Created on Sun Jun 28 18:31:23 2015

@author: hoseung
"""

def nout2lbt(nout, nout_fi=187):
    """
      A very simple function assuming delta a = 0.005,
      and nout_fi = 187 by default.
    """
    import astropy.cosmology as ac
    aexp = 1 - (nout_fi - nout)*0.005
    
    return ac.WMAP7.lookback_time(1/aexp -1).value



def time2gyr(times, z_now=None, info=None):
    """
    
    returns the age of "universe" at the given time.

    """
    import pyfits
    import numpy as np
    
#    repodir = '/home/hoseung/Copy/pyclusterevol/repo/'
    import inspect
    import os
    import utils
    repodir = os.path.dirname(inspect.getfile(utils)).split('/utils')[0] + '/repo/'

    if z_now is None and info is not None:
        z_now = 1/info.aexp-1
    z_now = max([z_now, 1e-10])
    sh0       = str(round(info.H0))
    som       = str(round(info.om*100))
    sol       = str(round(info.ol*100))

    tablefile  = repodir+'Table_taz_H'+sh0+'_Om'+som+'_Ol'+sol+'.fits'

    hdu = pyfits.open(tablefile)
    ttable = hdu[1].data

    tu       = ttable['t_unit'][0]
    tlb      = ttable['t_lback'][0]
    zred     = ttable['z'][0]
    aexp     = ttable['aexp'][0]

    t_lback_now = np.interp(z_now, zred, tlb)  # interpolation

    fd = np.where(times < min(tu))[0]
    if len(fd) > 0:
        ctime2 = times
        ctime2[fd] = min(tu)
        t_lback_in  = np.interp(ctime2, tu, tlb)
    else:
        t_lback_in  = np.interp(times, tu, tlb)

    return t_lback_in - t_lback_now


"""    
def friedman(O_mat_0, O_vac_0, O_k_0, alpha, ntable):
    if ((O_mat_0 + O_vac_0 + O_k_0) != 1.0):
        print('Error: wrong cosmological constants')
        print('O_mat_0,O_vac_0,O_k_0= ',O_mat_0, O_vac_0, O_k_0)
        print('The sum should be equal to 1.0, but ')
        print('O_mat_0+O_vac_0+O_k_0= ',O_mat_0 + O_vac_0 + O_k_0)
   
    axp_tau  = 1.0
    axp_t    = 1.0
    tau      = 0.0
    tau_prec = 0.0
    t        = 0.0
    t_prec   = 0.0
    ncount   = 0 
    
    print('tau, tau_out(last)', tau, tau_out(ntable-1))
    while (tau >= tau_out(ntable-1)):
        dtau        = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
        axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.0
        axp_tau     = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
        tau         = tau - dtau
    
        dt          = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
        axp_t_pre   = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.0
        axp_t       = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
        t           = t - dt
    
    if (tau <= tau_out[ncount] & tau_prec > tau_out[ncount]):
        t_out[ncount] = -t_prec
        zred_out[ncount] = 1.0/axp_t-1.0
        if (ncount < ntable-1):
            ncount = ncount + 1 
        print( ncount, ntable)

    tau_prec    = tau 
    t_prec      = t 
       
    print(' Age of the Universe (in unit of 1/H0)=',-t,-tau)
"""    
