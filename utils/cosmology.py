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



class Timeconvert():
    def __init__(self, info):
        import numpy as np
        from general import defaults
        from astropy.io import fits
        #import inspect
        #import os
        #import utils
        dfl = defaults.Default()
        self.repodir = dfl.dir_repo
        self.info = info

        sh0       = str(round(info.H0))
        som       = str(round(info.om*100))
        sol       = str(round(info.ol*100))

        tablefile  = self.repodir+'Table_taz_H'+sh0+'_Om'+som+'_Ol'+sol+'.fits'

        hdu = fits.open(tablefile)
        ttable = hdu[1].data

        self.tu       = ttable['t_unit'][0]
        self.tlb      = ttable['t_lback'][0]
        self.zred     = ttable['z'][0]
        self.aexp     = ttable['aexp'][0]

    def time2gyr(self, times, z_now=None):
        """
        
        returns the age of "universe" at the given time.

        """
        import numpy as np
        if z_now is None:
            z_now = 1/self.info.aexp-1

        z_now = max([z_now,1e-10])

        t_lback_now = np.interp(z_now, self.zred, self.tlb)  # interpolation

        fd = np.where(times < min(self.tu))[0]
        if len(fd) > 0:
            ctime2 = times
            ctime2[fd] = min(self.tu)
            t_lback_in  = np.interp(ctime2, self.tu, self.tlb)
        else:
            t_lback_in  = np.interp(times, self.tu, self.tlb)

        return t_lback_in - t_lback_now

