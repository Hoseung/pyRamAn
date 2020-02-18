# -*- coding: utf-8 -*-
"""
cosmology utils.

... use astropy.cosmology. that is a full furnished util.

Created on Sun Jun 28 18:31:23 2015

@author: hoseung
"""
from ..general import defaults
dfl = defaults.Default()
dir_repo = dfl.dir_repo

import numpy as np
class Timeconvert():
    def __init__(self, info=None, H0=None, om=None, ol=None, zred_now=None):
        from astropy.io import fits
        self.repodir = dir_repo
        self.info = info
        if info is not None:
            sh0       = str(round(info.H0))
            som       = str(round(info.om*100))
            sol       = str(round(info.ol*100))
            zred_now = info.zred
        else:
            sh0       = str(round(H0))
            som       = str(round(om*100))
            sol       = str(round(ol*100))
            zred_now  = zred_now

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
        self.t_lback_now = np.interp(zred_now, self.zred[::-1], self.tlb[::-1])  # interpolation

    def time2gyr(self, times, z_now=None):
        """
        returns the age of "universe" at the given time.
        """
        if z_now is not None:
            #z_now = max([z_now,1e-10])
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
            #z_now = max([z_now,1e-10])
            t_lback_now = np.interp(z_now, self.zred[::-1], self.tlb[::-1])
        else:
            t_lback_now = self.t_lback_now
        #
        t_lback_in  = np.interp(zreds, self.zred[::-1], self.tlb[::-1])
        return t_lback_in - t_lback_now
