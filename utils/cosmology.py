# -*- coding: utf-8 -*-
"""
cosmology utils.

... use astropy.cosmology. that is a full furnished util. 

Created on Sun Jun 28 18:31:23 2015

@author: hoseung
"""

import numpy as np
def nout2lbt(nout, nout_fi=187):
    """
      A very simple function assuming delta a = 0.005,
      and nout_fi = 187 by default.
    """
    import astropy.cosmology as ac
    aexp = 1 - (nout_fi - nout)*0.005
    
    return ac.WMAP7.lookback_time(1/aexp -1).value


class Timeconvert():
    def __init__(self, info):
        from general import defaults
        from astropy.io import fits
        dfl = defaults.Default()
        self.repodir = dfl.dir_repo
        self.info = info

        sh0       = str(round(info.H0))
        som       = str(round(info.om*100))
        sol       = str(round(info.ol*100))

        tablefile  = self.repodir+'Table_taz_H'+sh0+'_Om'+som+'_Ol'+sol+'.fits'

        hdu = fits.open(tablefile)
        ttable = hdu[1].data

        # Sort so that self.tu is in increasing order
        # Because converting stellar conformal times to lookback time 
        # is the main use case.
        # However, this sorting makes zred be a decreasing function.
        # So is needed the [::-1] indexing.
        isort=np.argsort(ttable['t_unit'][0])
        self.zred     = ttable['z'][0][isort]
        self.tu       = ttable['t_unit'][0][isort]
        self.tlb      = ttable['t_lback'][0][isort]
        self.aexp     = ttable['aexp'][0][isort]
        self.t_lback_now = np.interp(info.zred, self.zred[::-1], self.tlb[::-1])  # interpolation

    def time2gyr(self, times, z_now=None):
        """
        returns the age of "universe" at the given time.
        """
        if z_now is not None:
            z_now = max([z_now,1e-10])
            t_lback_now = np.interp(z_now, self.zred[::-1], self.tlb[::-1])
        else:
            t_lback_now = self.t_lback_now

        fd = np.where(times < min(self.tu))[0]
        if len(fd) > 0:
            ctime2 = times
            ctime2[fd] = min(self.tu)
            t_lback_in  = np.interp(ctime2, self.tu, self.tlb)
        else:
            t_lback_in  = np.interp(times, self.tu, self.tlb)

        return t_lback_in - t_lback_now

    def zred2gyr(self, zreds, z_now=None):
        if z_now is not None:
            z_now = max([z_now,1e-10])
            t_lback_now = np.interp(z_now, self.zred[::-1], self.tlb[::-1])
        else:
            t_lback_now = self.t_lback_now

        # zreds[zreds < 0] = 0 
        # 
        t_lback_in  = np.interp(zreds, self.zred[::-1], self.tlb[::-1])
        return t_lback_in - t_lback_now
        
        






