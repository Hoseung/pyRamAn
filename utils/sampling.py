# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 14:51:42 2015

Routines used in sampling tasks


@author: hoseung
"""
import numpy as np
class Region():
    def __init__(self, halo=None, **region_arguments):
        if halo is not None:
            self.region_from_halo(halo)
        else:
            self._xc = 0.5
            self._yc = 0.5
            self._zc = 0.5
            self._xr = [0,1]
            self._yr = [0,1]
            self._zr = [0,1]
            self._centers = [0.5, 0.5, 0.5]
            self._radius = 0.5
            self._ranges = [[0,1]]*3
            for key in region_arguments:
                try:
                    setattr(self, key, region_arguments[key])
                except:
                    continue
            #self.set_region(**region)

    @property
    def xc(self):
        return self._xc
    @xc.setter
    def xc(self, val):
        self._xc = val
        self._centers[0] = self._xc
        self._xr = [self._xc-self._radius,self._xc+self._radius]
        self._ranges[0] = self._xr

    @property
    def yc(self):
        return self._yc
    @yc.setter
    def yc(self, val):
        self._yc = val
        self._centers[1] = self._yc
        self._yr = [self._yc-self._radius,self._yc+self._radius]
        self._ranges[1] = self._yr

    @property
    def zc(self):
        return self._zc
    @zc.setter
    def zc(self, val):
        self._zc = val
        self._centers[2] = self._zc
        self._zr = [self._zc-self._radius,self._zc+self._radius]
        self._ranges[2] = self._zr

    @property
    def centers(self):
        return self._centers
    @centers.setter
    def centers(self, value):
        self._centers = value
        # update other values
        self._xc, self._yc, self._zc = self._centers
        self._xr = [self._xc-self._radius,self._xc+self._radius]
        self._yr = [self._yc-self._radius,self._yc+self._radius]
        self._zr = [self._zc-self._radius,self._zc+self._radius]
        self._ranges = [self._xr, self._yr, self._zr]

    @property
    def ranges(self):
        return self._ranges
    @ranges.setter
    def ranges(self, value):
        self._ranges = value
        self._xr, self._yr, self._zr = self._ranges
        self._xc, self._yc, self._zc = [0.5*sum(x) for x in self._ranges]
        self._centers = [self._xc, self._yc, self._zc]
        self._radius = max([self._xr[1] - self._xr[0],
                            self._yr[1] - self._yr[0],
                            self._zr[1] - self._zr[0]])

    @property
    def xr(self):
        return self._xr
    @xr.setter
    def xr(self, val):
        self._xr = val
        self._xc = 0.5*sum(self._xr)
        self._centers[0] = self._xc
        self._radius = max((self._radius, self._xr[1]-self._xr[0]))

    @property
    def yr(self):
        return self._yr
    @yr.setter
    def yr(self, val):
        self._yr = val
        self._yc = 0.5*sum(self._yr)
        self._centers[1] = self._yc
        self._radius = max((self._radius, self._yr[1]-self._yr[0]))

    @property
    def zr(self):
        return self._zr
    @zr.setter
    def zr(self, val):
        self._zr = val
        self._zc = 0.5*sum(self._zr)
        self._centers[2] = self._zc
        self._radius = max((self._radius, self._zr[1]-self._zr[0]))


    @property
    def radius(self):
        return self._radius
    @radius.setter
    def radius(self, val):
        self._radius = val
        self._xr = [self._xc - self._radius, self._xc + self._radius]
        self._yr = [self._yc - self._radius, self._yc + self._radius]
        self._zr = [self._zc - self._radius, self._zc + self._radius]
        self._ranges = [self._xr, self._yr, self._zr]


    def region_from_halo(self, halo, ind=None, rscale = 1.0):
        """
        Set a region as a halo.

        If multiple halos are given, define an area encompassing all halos.
        """
        if ind is None:
            try:
                self.centers = [halo['x'],
                                halo['y'],
                                halo['z']]
                self.radius = halo['rvir'] * rscale
            except:
                self.set_region_multi(xc = halo.data['x'],
                                      yc = halo.data['y'],
                                      zc = halo.data['z'],
                                      radius = halo.data['rvir'] * rscale)
        else:
            if len(ind) > 1:
                self.set_region_multi(xc = halo.data['x'][ind],
                                      yc = halo.data['y'][ind],
                                      zc = halo.data['z'][ind],
                                      radius = halo.data['rvir'][ind] * rscale)
            else:
                self.centers = [halo.data['x'][ind],
                                halo.data['y'][ind],
                                halo.data['z'][ind]]
                self.radius = halo.data['rvir'][ind] * rscale

    def distance_to(self, xc, xx):
        return np.sqrt([(xc[0] - xx[0])**2 + (xc[1] - xx[1])**2 + (xc[2] - xx[2])**2])[0]

    def extract_halos_within(self, halos, i_center, scale=1.0, Mcut=1e5):
        '''
        Return indices of halos within SCALE * Rvir of the central halo.

        def extract_halos_within(halos, ind_center, scale=1.0)
        halos : halo finder output (single snapshot)
        ind_center : index of the central halo
        scale : multiplying factor to the Rvir of the central halo

        NOTE
        ----
        Additional halo selection criteria are likely to be added onward.
        Therefore it is better to return mask array (True, True, ..., False, True)
        than a subset of indices.

        Another "extract_halos_within" is in lambda_mp_xxx.py series.
        That function accepts an explicit radius cut to work with galaxies. (rvir of galaxy is much smaller than rvir of halo.)

        To do
        -----
        TreeMaker fields are : ['p'][0:3],
        whereas HaloMaker gives ['x'], ['y'], ['z'].
        standardize it.
        '''
        xc = halos['x'][i_center]
        yc = halos['y'][i_center]
        zc = halos['z'][i_center]
        rvir= halos['rvir'][i_center]

        xx = halos['x']
        yy = halos['y']
        zz = halos['z']
        #m = np.array(halos['m'])
        dd = distance_to([xc,yc,zc], [xx,yy,zz])

        return (dd < (rvir * scale)) * (np.array(halos['m']) > Mcut)


    def set_region_multi(self, xc=None, yc=None, zc=None, radius=None):
        """
        Return region class.

        Note
        ----
        Accepts single / multiple inputs.
        If scalar variables are given, set_region is directly called.
        """
        if isinstance(xc, float):
            self.centers = [xc, yc, zc]
            self.radius = radius)
        else:
            ranges = self.sum_region_multi(xc=xc, yc=yc, zc=zc, radius=radius)
            rr = max([np.ptp(ranges[0:2]),
                      np.ptp(ranges[2:4]),
                      np.ptp(ranges[4:6])]) * 0.5

            self.centers = [0.5 * sum(ranges[0:2]),
                                   0.5 * sum(ranges[2:4]),
                                   0.5 * sum(ranges[4:6])]
            self.radius = rr


    def sum_region_multi(self, xc=[], yc=[], zc=[], radius=[]):
        """
        Returns minimum and maximum of given ranges.

        Note
        ----
        """

        xmi=[]
        xma=[]
        ymi=[]
        yma=[]
        zmi=[]
        zma=[]
        for i in range(len(xc)):
            xmi.append(xc[i] - radius[i])
            xma.append(xc[i] + radius[i])
            ymi.append(yc[i] - radius[i])
            yma.append(yc[i] + radius[i])
            zmi.append(zc[i] - radius[i])
            zma.append(zc[i] + radius[i])

        return (min(xmi), max(xma), min(ymi), max(yma), min(zmi), max(zma))
