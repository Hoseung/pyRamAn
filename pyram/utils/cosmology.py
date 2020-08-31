# -*- coding: utf-8 -*-
"""
cosmology utils.

... use astropy.cosmology. that is a full furnished util.

Created on Sun Jun 28 18:31:23 2015

@author: hoseung
"""
from ..general import defaults
#dfl = defaults.Default()
#dir_repo = dfl.dir_repo
from scipy.integrate import cumtrapz
from numpy.core.records import fromarrays as fromarrays
import numpy as np
from ..constants import km, Mpc, ly, Gyr

class Timeconvert():
    def __init__(self, info=None, H0=None, Om=None, Ol=None, zred_now=None, aexp_now=None, n=5000, verbose=True):
        if zred_now and aexp_now:
            assert zred_now == 1-1/aexp_now, f"You gave both redshift ({zred_now:.2f}) and exapnsion factor ({aexp_now:.2f}) and they differ."
        elif zred_now:
            aexp_now = 1/(1+zred_now)
        elif aexp_now:
            zred_now = 1/aexp_now -1
        elif info:
            zred_now = info.zred
            aexp_now = info.aexp
        else:
            zred_now = 0
            aexp_now = 1
            print("warning... Defaulting zred to 0, aexp to 1")

        self.zred_now = max([0,zred_now])
        self.aexp_now = min([1,aexp_now])

        if H0 is None: H0 = info.H0
        if Om is None: Om = info.om
        if Ol is None: Ol = info.ol
        # Integrate manually because astropy cosmology calculation is too slow...
        aarr = np.linspace(0, 1, n)[1:] ** 2
        aarr_st = (aarr[:-1] + aarr[1:])/2
        duda = 1. / (aarr_st ** 3 * np.sqrt(Om * aarr_st ** -3 + Ol))
        dtda = 1. / (H0 * km * Gyr / Mpc * aarr_st * np.sqrt(Om * aarr_st ** -3 + Ol))
        aarr = aarr[1:]

        uarr = cumtrapz(duda[::-1], aarr[::-1], initial=0)[::-1]
        tarr = cumtrapz(dtda, aarr, initial=0)
        self.cosmo_table = fromarrays([1/aarr - 1, aarr, tarr, uarr], dtype=[('zred', 'f8'), ('aexp', 'f8'), ('age', 'f8'), ('tu', 'f8')])
        self.age = np.interp(aexp_now, self.cosmo_table['aexp'], self.cosmo_table['age'])
        if verbose:
            print('Age of the universe (now/z=0): %.3f / %.3f Gyr, z = %.5f' % (self.age, self.cosmo_table['tu'][-1], zred_now))

    def time2gyr(self, times, zred_now=None, aexp_now=None):
        """
        Converts code unit times to lookback time from zred_now.
        age_of_star = Timeconvert.time2gyr(star["time"], zred_now=info.zred)

        NOTE
        ----
        uage : age of the Universe

        """
        if not zred_now and not aexp_now:
            zred_now = self.zred_now
        else:
            zred_now = self.zred_now
            aexp_now = self.aexp_now
        
        fd = np.where(times < min(self.cosmo_table.tu))[0]
        if len(fd) > 0:
            ctime2 = times
            ctime2[fd] = min(self.cosmo_table.tu)
            uage  = np.interp(ctime2, self.cosmo_table.tu, self.cosmo_table.age)
        else:
            uage  = np.interp(times, self.cosmo_table.tu, self.cosmo_table.age)

        return uage# - np.interp(zred_now, self.cosmo_table.zred, self.cosmo_table.age)

    def zred2gyr(self, zreds, zred_now=None, aexp_now=None):
        if not zred_now and not aexp_now:
            zred_now = self.zred_now
        else:
            zred_now = self.zred_now
            aexp_now = self.aexp_now

        fd = np.where(zreds < min(self.cosmo_table.zred))[0]
        if len(fd) > 0:
            ctime2 = zreds
            ctime2[fd] = min(self.cosmo_table.zred)
            uage  = np.interp(ctime2, self.cosmo_table.zred, self.cosmo_table.age)
        else:
            uage  = np.interp(zreds, self.cosmo_table.zred, self.cosmo_table.age)
        
        return uage - np.interp(zred_now, self.cosmo_table.zred, self.cosmo_table.age)
