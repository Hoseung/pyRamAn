# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 23:45:55 2015
@author: hoseung
"""
# Galaxy object
#class Galaxy(halo.HaloMeta):
"""
Inheries halo.HaloMeta.
HaloMeta does not contain positional and kinematic information of the halo.
This is good because recalculated center of galaxy is
likely to deviate from the center of the halo.


ATTENTION!!

with RdYlDu color map, empty cells in the velocity map are yellow.
cells with no mass should be NaN, not 0.
Fix it

"""

import numpy as np
from load.info import get_minimal_info
#import  matplotlib.pyplot as plt

def print_large_number(q):
    if isinstance(q, (int, np.int, np.int16, np.int32, np.int64)):
        return("{:d}".format(q))
    elif q > 1e4:
        return("{:.3e}".format(q))
    else:
        return("{:.2f}".format(q))


def convert_catalog(gcat, pboxsize, unit="code"):
    """
    Default unit = code unit.
    (tree.halomodule converts units to code units on reading gcats)
    """
    #gcat['x'] = (gcat['x'] -0.5) * pboxsize
    #gcat['y'] = (gcat['y'] -0.5) * pboxsize
    #gcat['z'] = (gcat['z'] -0.5) * pboxsize
    gcat['r'] = gcat['r'] * pboxsize  #* info.cboxsize * 1e3
    gcat['rvir'] = gcat['rvir'] * pboxsize #* info.cboxsize * 1e3


class Meta():
    def __init__(self):
        self.reff = 0.0
        self.mstar = 0.0
        self.nstar = 0
        self.mgas = 0.0
        #self.lambda_arr = None
        self.lambda_r = 0.0
        #self.lambda_arr2 = None # reoriented values
        #self.lambda_r2 = 0.0 # reoriented values
        self.vmap = None
        self.pt=[]
        self.pq=[]
        self.mrj=0.0
        self.d2t=0.0
        self.nvec=None
        self.lvec=None
        self.debug=False

    # System position
    @property
    def xc(self):
        return self._xc
    @xc.setter
    def xc(self, val):
        self._xc = val

    @property
    def yc(self):
        return self._yc
    @yc.setter
    def yc(self, val):
        self._yc = val

    @property
    def zc(self):
        return self._zc
    @zc.setter
    def zc(self, val):
        self._zc = val

    @property
    def pos(self):
        return (self._xc, self._yc, self._zc)
    @pos.setter
    def pos(self, val):
        self._xc, self._yc, self._zc  = val

    # System Velcity
    @property
    def vxc(self):
        return self._vxc
    @vxc.setter
    def vxc(self, val):
        self._vxc = val

    @property
    def vyc(self):
        return self._vyc
    @vyc.setter
    def vyc(self, val):
        self._vyc = val

    @property
    def vzc(self):
        return self._vzc
    @vzc.setter
    def vzc(self, val):
        self._vzc = val

    @property
    def vel(self):
        return (self._vxc, self._vyc, self._vzc)
    @vel.setter
    def vel(self, val):
        self._vxc, self._vyc, self._vzc  = val

    def show_summary(self):
        """
            prints a summary of the galaxy meta data.
        """

        for name in sorted(self.__dict__):
            q = self.__dict__[name]
            if q is not None:
                if isinstance(q, str):
                    print(name, ":", q)
                else:
                    try:
                        len(q) # if sequence
                        print(name, ":", q)
                    except:
                        if np.isfinite(q): # string..
                            print(name, ":", print_large_number(q))
                        else:
                            print(name, ":", q)
            else:
                print(name, ":", q)


class Galaxy():
    """
        Maybe.. this class is growing too heavy...?
        If I have to deal with millions of these...

        Case 1. Gcat only
        Case 2. Gcat and hcat
        Case 3. Hcat Only
        Case 4. Gal file Only
    """
    def __init__(self, info=None,
                 gcat=None, hcat=None, convert_cat=True):
        self.meta = Meta()
        self.set_info(info)
        if gcat is not None: self.set_gcat(gcat.copy(), convert_cat)
        if hcat is not None: self.set_hcat(hcat.copy(), convert_cat)
        self._has_star=False
        self._has_dm=False
        self._has_cell=False
        #self.meta.id=-1 # why did I do this...??

    def set_info(self, info):
        self.info = get_minimal_info(info)

    def set_substructure(self, gcat):

        self.gcat

    def set_hcat(self, halo, convert):
        """
            Copy so not to be affected by the original data being modified outside.
            center of galaxy may be updated by more sophisticated method.
        """
        is_tree=False
        if "pos" in halo.dtype.fields:
            is_tree = True
            tt = halo.copy()
            self.hcat = dict( x=tt["pos"][0],#*self.info.pboxsize,
            	      	  	y=tt["pos"][1],#*self.info.pboxsize,
            	      	 	z=tt["pos"][2],#*self.info.pboxsize,
            	     	  	vx=tt["vel"][0],
            	 	   	  	vy=tt["vel"][1],
            	  	        vz=tt["vel"][2],
                            r=tt["rvir"],
                            m=tt["m"]*1e11,
            	      	    id=tt["id"],
                            mvir=tt["mvir"],
                            rvir=tt["rvir"],
                            idx=tt["idx"])

            if not is_tree and convert:
                convert_catalog(self.hcat, self.info.pboxsize)

    def set_gcat(self, gcat, convert):
        """
            Copy so not to be affected by the original data being modified outside.
            center of galaxy may be updated by more sophisticated method.
        """
        is_tree=False
        # If a tree is given
        if "xp" in gcat.dtype.fields:
            print(gcat["xp"])
            if len(gcat["xp"]) == 3:
        	   	is_tree = True
        	   	tt = gcat.copy()
        	   	gcat = dict( x=tt["xp"][0],#*self.info.pboxsize,
              	      	   	  	y=tt["xp"][1],#*self.info.pboxsize,
              	      	   	 	z=tt["xp"][2],#*self.info.pboxsize,
              	      	   	  	vx=tt["vp"][0],
              	      	   	  	vy=tt["vp"][1],
              	      	        vz=tt["vp"][2],
                                r=tt["rvir"],
              	      	   	  	m=tt["m"]*1e11,
              	      	        id=tt["id"],
                                mvir=tt["mvir"],
                                rvir=tt["rvir"],
 	   	                        idx=tt["idx"])

        self.gcat = gcat

        if not is_tree and convert:
            convert_catalog(self.gcat, self.info.pboxsize)
        if gcat is not None:
            self.meta.id = int(gcat['id'])
            self.meta.xc = float(gcat["x"])
            self.meta.yc = float(gcat["y"])
            self.meta.zc = float(gcat["z"])


    def physical_print(self, x):
        print([i * self.info.pboxsize * 100 for i in list(x)])


    def _convert_unit_meta(self, unit_conversion):
        """
            convert meta data into convenient unit.
        """
        if unit_conversion == "code":
            self.pos *= self.info.pboxsize*1000
        elif unit_conversion == "GM":
            self.pos *= self.info.pboxsize
        else:
            print("[galaxy._convert_unit_meta] Unknown unit_conversion option")


    def _convert_unit(self, pop, unit_conversion):
        if unit_conversion is None:
            print("No target unit is specified")
            return

        data = getattr(self, pop)
        vxc_before = self.meta.vxc
        self.meta.vel = np.median(data["vel"], axis=0) * self.info.kms
        # Unit conversion
        if unit_conversion == "code":
            if "x" in data.dtype.names:
                data["pos"] = (data["pos"] - self.meta.pos) * self.info.pboxsize*1000
            if 'm' in data.dtype.names:
                data['m'] = data['m'] * self.info.msun
            if "vx" in data.dtype.names:
                print("vcx after", self.meta.vxc)
                print("vcx diff", self.meta.vxc - vxc_before)
                data['vel'] = data['vel'] * self.info.kms - self.meta.vel
        elif unit_conversion == "GM":
            if "x" in data.dtype.names:
                data['pos'] = (data['pos'] - self.meta.pos) * self.info.pboxsize*1e3
            if 'm' in data.dtype.names:
                data['m'] = data['m'] * 1e11 # in Msun.
            if "vx" in data.dtype.names:
                data['vel'] -= self.meta.vel

    def _add_dm(self, dm, idm):
        ndm_tot = len(idm)
        self.dm = np.recarray(ndm_tot, dtype=dm.dtype)
        if 'id' in dm.dtype.names:
            self.dm['id'] = dm['id'][idm]
        if 'm' in dm.dtype.names:
            if unit_conversion == "code":
                self.dm['m'] = dm['m'][idm] * self.info.msun
            elif unit_conversion == "GM":
                self.dm['m'] = dm['m'][idm] * 1e11

        if 'x' in dm.dtype.names: self.dm['pos'] = dm['pos'][idm]
        if 'vel' in dm.dtype.names: self.dm['vel'] = dm['vel'][idm]

        if verbose: print("DM data stored")


    def reff_main_gal(self, x, y, z):
        """
            Check for additional galaxy inside the region.
            If there is another galaxy, Rgal is limited to the
            Finally returns temporary galaxy radius.

            #. peaks of components are far enough.

            #. minor component is larger than `10%' of the major component (number of stars)
            * second peak is also a considerable galaxy, large enough to be considered as
            a minor merging counter part. (or minor interaction counter part.)

            #. minimum value is significantly smaller than the second peak.
            * two galaxise are not closely connected.
            * minor components are not background from a larger structure.

            Because smooth background should not have significant peak.
            Or, a background with significant peak would be the host galaxy of this galaxy.
            If that is the case, this galaxy should not be identified as a 'main' galaxy,
            and accounting for embeded satellites is another thing.
        """

        def _find_maxima(a, exclude_end=False):
            """
                Return indices of points that has larger values than neighbouring points
            """
            result = np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]
            if exclude_end:
                result[0] = False
                result[-1] = False
            return result

        def _find_minima(a, exclude_end=False):
            """
                Return indices of points that has smaller values than neighbouring points
            """
            result = np.r_[True, a[1:] < a[:-1]] & np.r_[a[:-1] < a[1:], True]
            if exclude_end:
                result[0] = False
                result[-1] = False
            return result

        def _smooth(y, box_pts):
            """
                Smooth array by convolution. (Center shift)

                Parameters
                ----------
                box_pts : int
                    size of convolution box
            """
            box = np.ones(box_pts)/box_pts
            y_smooth = np.convolve(y, box, mode='same')
            return y_smooth

        d = np.sqrt(x**2 + y**2 + z**2)
        nbins = 25
        n, bins = np.histogram(d, bins=nbins, normed=1)

        imax = np.arange(len(n))[_find_maxima(_smooth(n,3), exclude_end=True)]
        # Because I don't allow the first bin be a maximum, the number of bins
        # must be large enough so that the main galaxy is divided by enough bins.
        # With too few bins, the first bin, covering most of the galaxy by itself,
        # may have the largest number of particles.

        imin = np.arange(len(n))[_find_minima(_smooth(n,3), exclude_end=True)]

        # If there is no second component, imin3 =[]
        if len(imin) == 0 or len(imax) == 1:
            print("No second component.")
            print("Good to go!, Calcualte R_gal with all the particles")
            return False #self.get_radius(xx, yy, zz, mm, self.radius_method)

        else:
            cnt_major = sum(n[0:imin[0]])
            cnt_minor = sum(n[imin[0]:imax[1]])
            f_mM = cnt_minor/cnt_major
            print("particle number fraction, #star_Minor/#star_Major:", f_mM)
            if f_mM > 0.1:
                print("There is another galaxy at {:.2f}kpc from the center.".format(bins[imax[1]]))
                print("Suggested R_gal < {:.2f}kpc".format(bins[imin[0]]))
                return bins[imin[0]]
            else:
                return False

    def get_radius(self, xx, yy, zz, mm, method):
        """
            Method 1: Radius = Rvir_halo
            Method 2: Radius = half mass radius of given stellar particles
        """

        import numpy as np
        if method == "simple":
            return self.gcat['rvir']
        elif method == "eff":
            # requires galaxy center to be determined.
            dd = xx**2 + yy**2 + zz**2
            dmax = max(dd)
            r_again = dmax
            m_annuli=[]
            nbin=8
            for i in range(nbin):
                m_annuli.append(sum(
                mm[(dd > i * dmax/nbin) * (dd < (i+1)*dmax/nbin)]))

            for i in range(nbin):
                if m_annuli[i] > m_annuli[i-1]:
                    r_again = dmax * (i+1)/nbin
                    break
            i_again = dd < r_again
            dsort = np.argsort(dd[i_again])
            m = mm[i_again]
            # half-mass radius, not half-particle radius.
            cumsum = np.cumsum(m[dsort])

            ihalf = np.searchsorted(cumsum, 0.5 * cumsum[-1])
            # If two galaxies are merging,
            # Chances are cumsum has a flat interval.
            # If there is such anormaly, measure Reff again by
            # shirinking the radius a bit.
            reff = np.sqrt(dd[i_again][dsort[ihalf]])
            # Rough profile

            return reff

    def get_center(self, x, y, z, m, method='default', tol=1e-7, niter=10):
        """
            Determines the center of mass of given particles.
            IDL version iteratively searches for mass weighted mean position of
            particles within gradually decreasing search radius.
            Let's do the same thing here.
        """

        import numpy as np
        r = 10. / 200 * 1000  # 10kpc/h
        # Stops if r < 10pc
        # or, if dx_3d < 5 pc
        msum = sum(m)
        xc_old = sum(x * m)/msum
        yc_old = sum(y * m)/msum
        zc_old = sum(z * m)/msum

        i = 0
        while i < niter and r > tol :

            msum = sum(m)
            xc = sum(x * m)/msum
            yc = sum(y * m)/msum
            zc = sum(z * m)/msum

            dx3 = np.sqrt((xc_old - xc)**2 + (yc_old - yc)**2 + (zc_old - zc)**2)
            if dx3 < 5e-8:
                break

            # shrink the search radius by half.
            ind = np.where( (x-xc)**2 + (y-yc)**2 + (z-zc)**2 < r**2)
            if len(ind) < 100:
                break
            x = x[ind]
            y = y[ind]
            z = z[ind]

            r = 0.2 * max([x.ptp(), y.ptp(), z.ptp()])

            xc_old = xc
            yc_old = yc
            zc_old = zc

            i +=1

        return xc, yc, zc


    def get_center2(self, x,y,z,m, tol=3, max_iter=10):
        dx = 1.0
        def _ind2pos(i_this, ptp, dims):
            x = i_this[0]/dims[0]*ptp
            y = i_this[1]/dims[1]*ptp
            z = i_this[2]/dims[2]*ptp

            return x,y,z

        def _normalize(x):
            """
                Normalize an array into 0 <= array < 1.0
                by dividing by the array.ptp() and then
                moving the element with maximum value inside.
            """

            return (x - x.min())/x.ptp()


        def _get_dpeak(x,y,z,m,npix=9):
            from utils import assign
            ptp = max([x.ptp(), y.ptp(), z.ptp()])

            xmin, ymin, zmin = x.min(),y.min(),z.min()
            x = (x - xmin)/ptp * npix
            y = (y - ymin)/ptp * npix
            z = (z - zmin)/ptp * npix

            field = assign.cic(m, x, npix,
                               y, npix, z,
                               npix, wraparound=True,
                               average=False, norm_integer=False)
            i_peak = np.argmax(field) #
            a = i_peak // (npix*npix)
            b = (i_peak - a*npix*npix)//npix
            c = i_peak - a*npix*npix - b*npix
            xnew, ynew, znew =_ind2pos([c,b,a], ptp, [npix]*3)
            xnew = xmin + xnew
            ynew = ymin + ynew
            znew = zmin + znew

            return xnew, ynew, znew, dx

        for i in range(max_iter):
            xc_new, yc_new, zc_new, dx = _get_dpeak(x,y,z,m,npix=5)
            if dx < tol:
                return xc_new, yc_new, zc_new
            ind = np.where( (x-xc_new)**2 + (y-yc_new)**2 + (z-zc_new)**2 < dx**2)[0]
            x,y,z,m = x[ind],y[ind],z[ind],m[ind]

        return xc_new, yc_new, zc_new


    def e_vdiff(self, iv, vxt, vyt, vzt):
        vdiff = np.sqrt(((vxt - vxt[iv])**2 +
                         (vyt - vyt[iv])**2 +
                         (vzt - vzt[iv])**2))
        return sum(-1./vdiff[vdiff != 0])


    def get_cov(self, multithreaded=False, center_only=False):
        if center_only:
            # It is especially important when making galaxies
            # based on GalaxyMaker output because GalaxyMaker output
            # tends to be irregular.
            # Galaxies can have long tails.
            # So you need to consider only particles
            # near the main body of the galaxy.
            #
            ind = np.where(self.star['x']**2
                          +self.star['y']**2
                          +self.star['z']**2 < self.meta.reff )

            vx = self.star['vx'][ind]
            vy = self.star['vy'][ind]
            vz = self.star['vz'][ind]
            mm = self.star['m'][ind]
        else:
            vx = self.star['vx']
            vy = self.star['vy']
            vz = self.star['vz']

        vhalx = self.gcat['vx']
        vhaly = self.gcat['vy']
        vhalz = self.gcat['vz']

        r_frac = 0.7
        vrange = max([np.std(vx),
                      np.std(vy),
                      np.std(vz)]) * 3

        vind = np.where((vx - vhalx)**2 +
                        (vy - vhaly)**2 +
                        (vz - vhalz)**2 < (r_frac * vrange)**2)[0]

        npart = len(vind)
        method="com"
        vxt = vhalx
        vyt = vhaly
        vzt = vhalz
        # Add adaptivity
        if method == "mostbound":
            while npart > 80000:
                r_frac *= 0.7
                vind = np.where(np.square(vx - vxt) +
                    np.square(vy - vyt) +
                    np.square(vz - vzt) < r_frac * vrange**2)[0]
                npart = len(vind)
            vind = vind[0:-1:2]
            if npart < 100:
                return False
# sparsely selected representative particles.
            vxt = vx[vind]
            vyt = vy[vind]
            vzt = vz[vind]
            e = np.zeros(npart, dtype=float)
            vdiff = np.zeros(npart, dtype=float)

            # this takes ~ 30s. with roughly the same number of particles.
            for iv in range(1, len(vxt)-2):
                vdiff=np.sqrt(np.square(vxt[0:-1:5] - vxt[iv]) +
                              np.square(vyt[0:-1:5] - vyt[iv]) +
                              np.square(vzt[0:-1:5] - vzt[iv]))
                e[iv] = sum(-1./vdiff[vdiff != 0])

# N parts that are most bound to <50000 'central' particle group.
            emin = np.argsort(e)[:min(5000, npart)]

            mm = self.star['m'][vind]
            self.meta.vxc = np.average(vxt[emin], weights=mm[emin]) # top 100
            self.meta.vyc = np.average(vyt[emin], weights=mm[emin])
            self.meta.vzc = np.average(vzt[emin], weights=mm[emin])
        elif method == 'com':
            self.meta.vxc, self.meta.vyc, self.meta.vzc = self.get_center(
                                vx, vy, vz, mm, tol=5, niter=8)

    def cal_mgas(self):
        #msun = 1.98892e33 # solar mass in gram.
        self.meta.mgas = sum(self.cell['rho'] * self.cell['dx']**3) * self.info.punit_m
        # [g/cm^3] * [cm^3] / [g/msun] = [msun]

    def cal_ssfr(self):
        if self.meta.sfr == 0:
            self.cal_sfr()
        self.meta.ssfr = self.meta.sfr / self.meta.mstar

    def dist_map(self, npix=40):
        nx = npix
        ny = npix

        dist_map=np.zeros((nx, ny), dtype=float)
        for i in range(nx):
            for j in range(ny):
                dist_map[i][j]= np.sqrt((0.5 + i - nx/2)**2 + (0.5 + j - ny/2)**2)
        return dist_map

    def weighted_std(self, values, weights):
        import math
        """
        Return the weighted average and standard deviation.

        values, weights -- Numpy ndarrays with the same shape.
        """

        average = np.average(values, weights=weights)
        variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise

        return math.sqrt(variance)


    def _pseudo_particles(self,xx,yy,mm,vz, sig=0.5, n_times=60):
            cov = [[sig,0],[0,sig]]
            dx, dy = np.random.multivariate_normal([0,0], cov, len(xx) * n_times).T
            xnew = np.repeat(xx, n_times) + dx
            ynew = np.repeat(yy, n_times) + dy
            mm = np.repeat(mm,n_times)/n_times
            vz = np.repeat(vz,n_times)
            return xnew, ynew, mm, vz


    def reorient(self, dest=[0., 0., 1], pop_nvec = ['star'],
                 pops=['star', 'dm', 'cell'], verbose=False, debug=False):
        """
        rotate particles/cells using rotation matrix.

        Parameters
        ----------
        dest : float array [3]
            direction vector the normal vector will point to.

        pop_nvec : a list of galaxy components ['star', 'dm', 'cell']
            target_nvec determines with respect to which the galaxy rotates.
            stellar particels / DM particles/ gas / all particles.
            Defaults to ['star'].
        verbose :
            verbosity.

        Returns
        -------
            Nothing

        Notes
        -----
            Use stellar particles by default.

        Examples
        --------

        See Also
        --------

        References
        ----------

        """
        if not hasattr(self.meta, 'nvec'):
            """
                fix : nvec w.r.t what?, and what is re-oriented?
            """
            if verbose:
                print("No normal vector, calculating for {} population".format(pop_nvec))
            self.cal_norm_vec(pop_nvec)

        if not hasattr(self, 'rotation_matrix'):
            if verbose:
                print("Galaxy.reorientation: No Rotation matrix. \n Calculating..")
            self.cal_rotation_matrix(dest = dest)

        if debug:
            print("---------------------")
            print(self.meta.nvec)
            print(self.rotation_matrix)
            print("---------------------")

        # New coordinates
        for target in pops:
            pop = getattr(self, target)
            RM = self.rotation_matrix
            pop["pos"] = np.matmul(pop["pos"], RM.T)# (pop['x'], pop['y'], pop['z'])))
            pop["vel"] = np.matmul(pop["vel"], RM.T)#, pop['vy'], pop['vz'])))

    def cal_norm_vec(self, pop_nvec=['star'], dest=[0., 0., 1], bound_percentile=50):
        """
        Determines normal vector (rotation axis) of the galaxy.

        Parameters
        ----------
        pop_nvec : ["name", "of", "species"]
            components for which normal vector is calculated.
        dest : float array [3]
            direction vector the normal vector will point to.
        bound_percentile :
            Top ?% bound particles are taken into calculation.

        Examples
        --------
        >>> self.cal_norm_vec(pop_nvec=['star'], dest=[0., 0., 1.], bound_percentile=50)

        Notes
        -----
        'galaxy' means stellar particles by default.
        Add self.nvec attribute and also returns nvec


        See Also
        --------
            Naab + 2014: normal vector from the most bound 50% particles.
        """

        nelements = 0

        for target in pop_nvec:
            pop = getattr(self, target)
            #nelements += len(pop['x'])

        if bound_percentile < 100:
            vdiff = np.sqrt(np.sum(np.square(self.star["vel"]), axis=1))

            i_bound_50 = np.argsort(vdiff)[:max([1, 0.01 * bound_percentile]) * len(vdiff)]

            #pop = getattr(self, target)
            #nn = len(pop['x'])
            pos = self.star['pos'][i_bound_50]
            vel = self.star["vel"][i_bound_50]

        else:
            """
                all species are taken into account.
            """
            nelements = sum(self.bound_ptcl)
            pos = np.zeros(nelements,3)
            vel = np.zeros(nelements,3)
            iskip = 0
            for target in pop_nvec:
                pop = getattr(self, target)
                nn = len(pop['x'])
                pos[iskip:iskip + nn] = pop['pos'][self.bound_ptcl]
                vel[iskip:iskip + nn] = pop['vel'][self.bound_ptcl]
                iskip += nn

        self.meta.lvec = np.sum(np.cross(pos,vel), axis=0)
        self.meta.nvec = self.meta.lvec/np.linalg.norm(self.meta.lvec)
        return self.meta.nvec

    def _get_rotation_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        http://stackoverflow.com/a/6802723
        https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
        """

        #import math
        axis = np.asarray(axis)
        theta = np.asarray(theta)
        axis = axis/np.sqrt(np.dot(axis, axis))
        a = np.cos(theta/2)
        b, c, d = -axis*np.sin(theta/2)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


    def cal_rotation_matrix(self, nvec=None, dest=[0., 0., 1]):
        """
        Determine Rotation Matrix by axis-angle rotation.

        parameters
        -------------
        nvec : normal vector.
            If not given, automatically calculated.
        dest : destination vector ( not axis of rotation!).
               +z direction ([0, 0, 1]) by default.
        """

        import numpy.linalg as lag
        import math
        # rotation axis
        dest = np.asarray(dest) # towards +z direction.
        if nvec is None:
            try:
                nvec = self.meta.nvec
            except:
                nvec = self.cal_norm_vec() # Calculate, add, and return .nvec attribute.

        if lag.norm(nvec) != 1.0:
            nvec = nvec / lag.norm(nvec)
        if lag.norm(dest) != 1.0:
            dest = dest / lag.norm(dest)

        print(nvec, dest)
        r_axis = np.cross(nvec, dest)
        angle = math.acos(np.dot(nvec, dest))

        self.rotation_matrix = self._get_rotation_matrix(r_axis, angle)
        return self.rotation_matrix

    def _follow_bp(self, bp):
        """
            Follow bound particle..

        """

    def cal_bound_ptcls(self, ptype='star', bound50=True):
        def gas_mass(cell, info):
            """
                return gas mass of cells.
                Keep in mind that the size of cells differ.
            """
            msun = 1.98892e33 # solar mass in gram.
            return (cell['rho'] * info.unit_d) * (cell['dx'] * info.unit_l)**3 / msun
        # what is bound particle?

        G = 6.67384e-11  # m^3 kg^-1 s^-2
        kpc_to_m = 3.08567758e19
        msun_to_kg = 1.9891e30 # kg

        m_s = self.star['m']# * info.msun
        m_g = gas_mass(self.cell, self.info)
        m_d = self.dm['m']# * info.msun

        r_s = np.sqrt(np.sum(np.square(self.star["pos"]), axis=1))
        r_g = np.sqrt(np.sum(np.square(self.cel["pos"]), axis=1))
        r_d = np.sqrt(np.sum(np.square(self.dm["pos"]), axis=1))

        m_all = np.concatenate((m_s, m_g, m_d))
        r_all = np.concatenate((r_s, r_g, r_d))

        i_sorted = np.argsort(r_all) # r_all[i_sorted] = 0, 0.1, 0.2, 0.3, ... 100
        m_enc = np.cumsum(m_all[i_sorted])

        i_star = np.searchsorted(r_all[i_sorted], r_s)

        v_bound = np.sqrt(2*G*m_enc[i_star] * msun_to_kg/ (r_s * kpc_to_m)) * 1e-3

        vdiff = np.sqrt(np.sum(np.square(self.star["vel"]), axis=1))

        self.bound_ptcl = vdiff < v_bound


    def bound_particles(self, n=100, ptype='star'):
        """
        save IDs of n most bound particles as galaxy.most_bound_particles.

        Parameters
        ----------
        n : int
            number of most bound particles to be kept.
        ptype : str "star" or "dm"
            Type of particle.

        NOTE
        ----
        The result has no use yet... I don't know why I have this...
        Requires center of mass and center of velocity already determined.
        """

        #vrange = 100
        vx = self.star['vx']
        vy = self.star['vy']
        vz = self.star['vz']

        vdiff = np.square(vx - self.meta.vxc) + np.square(vy - self.meta.vyc) + np.square(vz - self.meta.vzc)
        npart = len(vdiff)

        if npart < 100:
            print("Too little particles within velocity window km/s")
            return False

        vdiff = np.zeros(npart, dtype=float)

        emin = np.argsort(vdiff)[:min(n, npart)]
        self.most_bound_particles = self.star['id'][emin]


    def cal_b2t(self, ptype='star', disk_criterion="Abadi2003",
                bound_only = True, return_ellip=False):
        """
        Measure bulge to total ratio of the galaxy.

        parameters
        ----------
        ptype : str
            type of particle to consider.
        disk_criterion : str
            discrimination criterion. (Defaults to "Abadi2003")
        bound_only : logical
            consider bound particles only
        return_ellip : logical
            if True, ellipticity of each particle is returned
        """
        def gas_mass(cell, info):
            """
                return mass of gas cells.
                Keep in mind that size of cells differ.
            """
            msun = 1.98892e33 # solar mass in gram.
            return (cell['rho'] * info.unit_d) * (cell['dx'] * info.unit_l)**3 / msun

        #if proj=="z":
        pos1, pos2, vel1, vel2 = 'x', 'y', 'vx', 'vy'


        # Radius

        # Calculating boundness requires total mass inside a radius.
        # -> DM, Cell are also needed.
        #part = getattr(self, ptype)
        m_d = self.dm['m']# * info.msun

        if hasattr(self, "cell"):
            m_g = gas_mass(self.cell, self.info)
            m_all = np.concatenate((self.star['m'], m_g, m_d))
        else:
            m_all = np.concatenate((self.star['m'], m_d))

        r_s = np.sqrt(np.sum(np.square(self.star["pos"]), axis=1))
        r_d = np.sqrt(np.sum(np.square(self.dm["pos"]), axis=1))

        if hasattr(self, "cell"):
            r_g = np.sqrt(np.sum(np.square(self.cell["pos"]), axis=1))
            r_all = np.concatenate((r_s, r_g, r_d))
        else:
            r_all = np.concatenate((r_s, r_d))

        i_sorted = np.argsort(r_all)
        m_enc = np.cumsum(m_all[i_sorted])

        # First nstar indices are stars.
        if bound_only:
            i_star = i_sorted[self.bound_ptcl]
            x = self.star[pos1][self.bound_ptcl]
            y = self.star[pos2][self.bound_ptcl]
            vx = self.star[vel1][self.bound_ptcl]
            vy = self.star[vel2][self.bound_ptcl]
            m = self.star['m'][self.bound_ptcl]# * info.msun
        else:
            i_star = i_sorted[0:len(r_s)]
            x = self.star[pos1]
            y = self.star[pos2]
            vx = self.star[vel1]
            vy = self.star[vel2]
            m = self.star['m']# * info.msun

        #boxtokpc = self.info.pboxsize * 1000
        G = 6.67384e-11  # m^3 kg^-1 s^-2
        kpc_to_m = 3.08567758e19
        msun_in_kg = 1.9891e30 # kg
        v_circ = np.sqrt(G * msun_in_kg * m_enc[i_star]/
                        (kpc_to_m * r_all[i_star])) * 1e-3 # m/s to in km/s

        j_circ = np.sqrt(np.square(x)
                    + np.square(y)) * v_circ # * boxtokpc omitted (later compensated.)
        # Finally, r in code unit, v in km/s
        j_phi = (x * vy - y * vx) # * boxtokpc omitted.
        ellipticity = j_phi / j_circ

        if disk_criterion == "Scannapieco2009":
            # < 0.8 from Scannapieco 2009
            disk = np.sum(m[ellipticity < 0.8])
        elif disk_criterion == "Abadi2003":
            #bulge = 2.0 * np.sum(self.star['m'][j_phi < 0])
            bulge = 2.0 * min([np.sum(m[j_phi < 0]),
                               np.sum(m[j_phi > 0])])
            disk = np.sum(m) - bulge

        self.meta.d2t = disk / self.meta.mstar # or sum(self.star['m'][ind])
        self.meta.b2t = bulge / self.meta.mstar
        if return_ellip:
            return ellipticity


    def cal_trivia(self):
        self.meta.mstar = sum(self.star['m'])
        self.meta.vcen = self.get_vcen()


    def save_gal_pickle(self):
        import pickle
        with open(str(self.meta.id).zfill(6) + 'gal.pickle', 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
