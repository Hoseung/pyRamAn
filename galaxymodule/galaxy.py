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


def convert_catalog(catalog, pboxsize, unit="code"):
    """
    Default unit = code unit.
    (tree.halomodule converts units to code units on reading catalogs)
    """

    catalog['x'] = (catalog['x'] -0.5) * pboxsize
    catalog['y'] = (catalog['y'] -0.5) * pboxsize
    catalog['z'] = (catalog['z'] -0.5) * pboxsize
    catalog['r'] = catalog['r'] * pboxsize  #* info.cboxsize * 1e3
    catalog['rvir'] = catalog['rvir'] * pboxsize #* info.cboxsize * 1e3


class Meta():
    def __init__(self):
        self.xc = 0.0
        self.yc = 0.0
        self.zc = 0.0
        self.vxc = 0.0
        self.vyc = 0.0
        self.vzc = 0.0
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
        #self.sfr=0.0
        #self.ssfr=0.0
        self.mrj=0.0
        self.d2t=0.0
        self.nvec=None
        self.lvec=None
        self.debug=False

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
    """
    def __init__(self, info=None, halo=None,
                 catalog=None, convert_cat=True):
        self.meta = Meta()
        #if info is not None:
        self.set_info(info)
        self.set_catalog(np.copy(catalog), convert_cat)
        self._has_star=False
        self._has_dm=False
        self._has_cell=False
        #self.meta.id=-1 # why did I do this...??

    def set_info(self, info):
        self.info = get_minimal_info(info)

    def set_catalog(self, catalog, convert):
        """
            Copy so not to be affected by the original data being modified outside.
            center of galaxy may be updated by more sophisticated method.
        """
        is_tree=False
        if "xp" in catalog.dtype.fields:
            print(catalog["xp"])
            if len(catalog["xp"]) == 3:
        	   	is_tree = True
        	   	tt = catalog.copy()
        	   	catalog = dict( x=tt["xp"][0],#*self.info.pboxsize,
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

        self.gcat = catalog

        if not is_tree and convert:
            convert_catalog(self.gcat, self.info.pboxsize)
        if catalog is not None:
            self.meta.id = int(catalog['id'])
            self.meta.xc = catalog["x"]
            self.meta.yc = catalog["y"]
            self.meta.zc = catalog["z"]

    def set_ptypes(self, pt=None):
        pass

    def physical_print(self, x):
        # If length,

        print([i * self.info.pboxsize * 100 for i in list(x)])


    def _convert_unit_meta(self, unit_conversion):
        """
            convert meta data into convenient unit.
        """
        if unit_conversion == "code":
            self.meta.xc *= self.info.pboxsize*1000
            self.meta.yc *= self.info.pboxsize*1000
            self.meta.zc *= self.info.pboxsize*1000
        elif unit_conversion == "GM":
            self.meta.xc *= self.info.pboxsize
            self.meta.yc *= self.info.pboxsize
            self.meta.zc *= self.info.pboxsize
        else:
            print("[galaxy._convert_unit_meta] Unknown unit_conversion option")


    def _convert_unit(self, pop, unit_conversion):
        if unit_conversion is None:
            print("No target unit is specified")
            return

        data = getattr(self, pop)
        print("vcx before", self.meta.vxc)
        vxc_before = self.meta.vxc
        self.meta.vxc = np.median(data["vx"]) * self.info.kms
        self.meta.vyc = np.median(data["vy"]) * self.info.kms
        self.meta.vzc = np.median(data["vz"]) * self.info.kms
        # Unit conversion
        if unit_conversion == "code":
            if "x" in data.dtype.names:
                data['x'] = (data['x'] - self.meta.xc) * self.info.pboxsize*1000
                data['y'] = (data['y'] - self.meta.yc) * self.info.pboxsize*1000
                data['z'] = (data['z'] - self.meta.zc) * self.info.pboxsize*1000
            if 'm' in data.dtype.names:
                data['m'] = data['m'] * self.info.msun
            if "vx" in data.dtype.names:
                print("vcx after", self.meta.vxc)
                print("vcx diff", self.meta.vxc - vxc_before)
                data['vx'] = data['vx'] * self.info.kms - self.meta.vxc
                data['vy'] = data['vy'] * self.info.kms - self.meta.vyc
                data['vz'] = data['vz'] * self.info.kms - self.meta.vzc
        elif unit_conversion == "GM":
            if "x" in data.dtype.names:
                data['x'] = (data['x'] - (self.meta.xc - 0.5) * self.info.pboxsize)*1e3
                data['y'] = (data['y'] - (self.meta.yc - 0.5) * self.info.pboxsize)*1e3
                data['z'] = (data['z'] - (self.meta.zc - 0.5) * self.info.pboxsize)*1e3
            if 'm' in data.dtype.names:
                data['m'] = data['m'] * 1e11 # in Msun.
            if "vx" in data.dtype.names:
                data['vx'] -= self.meta.vxc
                data['vy'] -= self.meta.vyc
                data['vz'] -= self.meta.vzc

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

        if 'x' in dm.dtype.names:
            self.dm['x'] = dm['x'][idm]
            self.dm['y'] = dm['y'][idm]
            self.dm['z'] = dm['z'][idm]

        if 'vel' in dm.dtype.names:
            self.dm['vx'] = dm['vx'][idm]
            self.dm['vy'] = dm['vy'][idm]
            self.dm['vz'] = dm['vz'][idm]

        if verbose: print("DM data stored")

    def _add_cell(self, cell, icell):
        ncell_tot = len(icell)
        dtype_cell = [('x', '<f8',),('y', '<f8',),('z', '<f8',), ('dx', '<f8'),
                      ('rho', '<f8'), ('vx', '<f8' ), ('vy', '<f8' ), ('vz', '<f8' ),
                      ('temp', '<f8'), ('metal', '<f8')]

        self.cell = np.recarray(ncell_tot, dtype=dtype_cell)
        self.cell['x'] = cell['x'][icell]
        self.cell['y'] = cell['y'][icell]
        self.cell['z'] = cell['z'][icell]
        self.cell['dx'] = cell['dx'][icell]
        self.cell['rho'] = cell['var0'][icell]
        self.cell['vx'] = cell['var1'][icell]
        self.cell['vy'] = cell['var2'][icell]
        self.cell['vz'] = cell['var3'][icell]
        self.cell['temp'] = cell['var4'][icell]
        self.cell['metal'] = cell['var5'][icell]


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

#            print(x.min(),x.max(),y.min(),y.max(),z.min(),z.max())
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
        #kpc_in_cm = 3.086e21
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
            By default, stellar particles by default.

        Examples
        --------

        See Also
        --------

        References
        ----------

        """
        if not hasattr(self, 'nvec'):
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
            new_x = np.matmul(RM, np.vstack((pop['x'], pop['y'], pop['z'])))
            new_v = np.matmul(RM, np.vstack((pop['vx'], pop['vy'], pop['vz'])))

            pop['x'] = new_x[0,:]
            pop['y'] = new_x[1,:]
            pop['z'] = new_x[2,:]

            pop['vx'] = new_v[0,:]
            pop['vy'] = new_v[1,:]
            pop['vz'] = new_v[2,:]

#            print("Error in reorient, couldn't update pos and vel")
#            continue

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
            vx = self.star['vx']
            vy = self.star['vy']
            vz = self.star['vz']
            vdiff = np.sqrt(np.square(vx) + np.square(vy) + np.square(vz))

            i_bound_50 = np.argsort(vdiff)[:max([1, 0.01 * bound_percentile]) * len(vdiff)]

            #pop = getattr(self, target)
            #nn = len(pop['x'])
            x = self.star['x'][i_bound_50]
            y = self.star['y'][i_bound_50]
            z = self.star['z'][i_bound_50]
            vx = vx[i_bound_50]
            vy = vy[i_bound_50]
            vz = vz[i_bound_50]

        else:
            """
                all species are taken into account.
            """
            nelements = sum(self.bound_ptcl)
            x = np.zeros(nelements)
            y = np.zeros(nelements)
            z = np.zeros(nelements)
            vx = np.zeros(nelements)
            vy = np.zeros(nelements)
            vz = np.zeros(nelements)
            iskip = 0
            for target in pop_nvec:
                pop = getattr(self, target)
                nn = len(pop['x'])
                x[iskip:iskip + nn] = pop['x'][self.bound_ptcl]
                y[iskip:iskip + nn] = pop['y'][self.bound_ptcl]
                z[iskip:iskip + nn] = pop['z'][self.bound_ptcl]
                vx[iskip:iskip + nn] = pop['vx'][self.bound_ptcl]
                vy[iskip:iskip + nn] = pop['vy'][self.bound_ptcl]
                vz[iskip:iskip + nn] = pop['vz'][self.bound_ptcl]
                iskip += nn

        lx = sum(y*vz - z*vy) # normalized; *m is omitted.
        ly = sum(z*vx - x*vz)
        lz = sum(x*vy - y*vx)
        self.meta.lvec = np.array([lx, ly, lz])
        self.meta.nvec = self.meta.lvec/np.sqrt(lx**2 + ly**2 + lz**2)
#        self.cal_rotation_matrix(dest=dest)
        return self.meta.nvec

    def _get_rotation_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        http://stackoverflow.com/a/6802723
        https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
        """

        import math
        axis = np.asarray(axis)
        theta = np.asarray(theta)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2)
        b, c, d = -axis*math.sin(theta/2)
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

        r_axis = np.cross(nvec, dest)
        angle = math.acos(np.dot(nvec, dest))

        RM = self._get_rotation_matrix(r_axis, angle)
        self.rotation_matrix = RM
        return RM

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

        r_s = np.sqrt(np.square(self.star['x'])
                    + np.square(self.star['y'])
                    + np.square(self.star['z']))
        r_g = np.sqrt(np.square(self.cell['x'])
                    + np.square(self.cell['y'])
                    + np.square(self.cell['z']))
        r_d = np.sqrt(np.square(self.dm['x'])
                    + np.square(self.dm['y'])
                    + np.square(self.dm['z']))

        m_all = np.concatenate((m_s, m_g, m_d))
        r_all = np.concatenate((r_s, r_g, r_d))

        i_sorted = np.argsort(r_all) # r_all[i_sorted] = 0, 0.1, 0.2, 0.3, ... 100
        m_enc = np.cumsum(m_all[i_sorted])

        i_star = np.searchsorted(r_all[i_sorted], r_s)

        v_bound = np.sqrt(2*G*m_enc[i_star] * msun_to_kg/ (r_s * kpc_to_m)) * 1e-3

        vx = self.star['vx']
        vy = self.star['vy']
        vz = self.star['vz']

        vdiff = np.sqrt(np.square(vx) + np.square(vy) + np.square(vz))

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
                bound_only = True, hist=False, proj="z"):
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
        hist : logical
            if True, ellipticity histogram is returned
        """
        def gas_mass(cell, info):
            """
                return mass of gas cells.
                Keep in mind that size of cells differ.
            """
            msun = 1.98892e33 # solar mass in gram.
            return (cell['rho'] * info.unit_d) * (cell['dx'] * info.unit_l)**3 / msun

        if proj=="x":
            pos1, pos2, vel1, vel2 = 'y', 'z', 'vy', 'vz'
        if proj=="y":
            pos1, pos2, vel1, vel2 = 'z', 'x', 'vz', 'vx'
        if proj=="z":
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


        r_s = np.sqrt(np.square(self.star['x'])
                    + np.square(self.star['y'])
                    + np.square(self.star['z']))

        if hasattr(self, "cell"):
            r_g = np.sqrt(np.square(self.cell['x'])
                    + np.square(self.cell['y'])
                    + np.square(self.cell['z']))
        r_d = np.sqrt(np.square(self.dm['x'])
                    + np.square(self.dm['y'])
                    + np.square(self.dm['z']))

        if hasattr(self, "cell"):
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
        if hist:
            return np.histogram(ellipticity, range=[-2,2], bins=20)



    def cal_trivia(self):
        self.meta.mstar = sum(self.star['m'])
        self.meta.vcen = self.get_vcen()

    def save_gal_pickle(self):
        import pickle
        with open(str(self.meta.id).zfill(6) + 'gal.pickle', 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)


def plot_gal(self, npix=200, fn_save=None, ioff=True,
             do9plots = False, i_lambda=0, **kwargs):
    """
        parameters
        ----------
        i_lambda: int
            index of gal.meta.lambda_result_list to be plotted as lambda profile.

    """

    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from draw import pp
    from matplotlib.colors import LogNorm

    if ioff:
        plt.ioff()

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)


    def _close_to(values, v_ref):
        """
            returns the index of the value which is close to the v_ref.
        """
        ind = 0
        good=values[ind]
        for i, vv in enumerate(values):
            if abs(vv - v_ref) < good:
                good = vv
                ind = i

        return ind

    rgal_to_reff = self.meta.Rgal_to_reff

    if do9plots:
        fig, axs = plt.subplots(3,3)
        fig.set_size_inches(12,9) # 9 plots
    else:
        fig, axs = plt.subplots(2,2)

    try:
        fig.suptitle("ID: {}    z: {:2f}".format(str(self.meta.id).zfill(5), self.info.zred))
    except:
        fig.suptitle("ID: {}    z: not available".format(str(self.meta.id).zfill(5)))

# Stellar particle density map
    ax = axs[0,0]
    """
        1st plot. - Stellar density map.

    """
    # if hist=True, use histogram instead of custom CIC.
    img = pp.part2den(self.star, self.info, npix=npix, hist=True)
    im = ax.imshow(img.data, cmap=plt.get_cmap('brg'),
                   norm=LogNorm(vmin=1e6))
#        _,_,_,im = ax.hist2d( self.star['x'], self.star['y'], norm=LogNorm(), bins=npix)
    fig.colorbar(im, ax=ax)

    ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
#   Todo:
#   If rgal_to_reff is not too large, say 10. xticks are too dense.
#   If 8, or 9 is given, use 4, if 10, use 5, and so on..
#   Utilize _close_to function.

    ax.set_xticks(np.linspace(0, npix, rgal_to_reff))
    xticks = np.linspace(-rgal_to_reff, rgal_to_reff, rgal_to_reff)
    ax.set_xticklabels(["{:.1f}".format(xx) for xx in xticks])

    ax.set_ylabel("[kpc]")
    ax.set_yticks(np.linspace(0,npix, 5))
    yticks = ["{:.2f}".format(y) \
                for y in np.linspace(-self.meta.rgal, self.meta.rgal, num=5)]
    ax.set_yticklabels(yticks)

# Lambda_r plot
    try:
        ll = self.meta.lambda_result_list[i_lambda]
        if ll is not None:
            ax = axs[0,1]
            ax.plot(ll) # ~ 1 * Reff
            ax.set_title(r"$\lambda _{R}$")
            ax.text(0.5 * len(ll), 0.8, "{:.2e}".format(self.meta.mstar) + r"$M_{\odot}$")
            ax.set_ylim(bottom=0, top=1)
            ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
 
            unit = len(ll)/self.meta.rscale_lambda # number of points per 1Reff
            xticks = [0, 0.5*unit]
            reffs=[0, 0.5]
            for i in range(int(self.meta.rscale_lambda)): # If lambda_r is calculated furthrer -
                xticks.append(unit * (i+1))
                reffs.append((i+1))
            ax.set_xticks(xticks)
            xticks_label = ["{:.1f}".format(rf) for rf in reffs]
            ax.set_xticklabels(xticks_label)
    except:
        pass

# sigma map
    if hasattr(self, "sigmap"):
        l_range = self.meta.reff * self.meta.rscale_lambda
        ax = axs[1,0]
        im = ax.imshow(self.sigmap, vmin=0, vmax = 200, cmap=plt.get_cmap('brg'))
        cb = plt.colorbar(im, ax=ax, label=r'$km s^{-1}$') # needs a mappable object.
        tick_locator = ticker.MaxNLocator(nbins=3)
        cb.locator = tick_locator
        cb.update_ticks()
        ax.set_title(r"$\sigma$ map")
        ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
        ax.set_xticks(np.linspace(0, len(self.mmap), num=3))
        ax.set_xticklabels([str(-1*self.meta.rscale_lambda),"0", str(self.meta.rscale_lambda)])
        ax.set_ylabel("[kpc]")
        ax.set_yticks(np.linspace(0, len(self.mmap), num=5))
        yticks = ["{:.2f}".format(y) \
                    for y in np.linspace(-l_range, l_range, num=5)]
        ax.set_yticklabels(yticks)
#        ax.locator_params(tight=True, nbins=5)

# velocity map
    if hasattr(self,"vmap") and self.vmap is not None:
        l_range = self.meta.reff * self.meta.rscale_lambda
        ax = axs[1,1]
        im = ax.imshow(self.vmap, vmin=-100, vmax = 100, cmap='RdBu')
#        ax.tick_params(
#            which='both',      # both major and minor ticks are affected
#            bottom='off',      # ticks along the bottom edge are off
#            top='off',         # ticks along the top edge are off
#            labelbottom='off')
        cb = plt.colorbar(im, ax=ax, label=r'$km s^{-1}$') # needs a mappable object.
        tick_locator = ticker.MaxNLocator(nbins=3)
        cb.locator = tick_locator
        cb.update_ticks()
        im.set_cmap('RdYlBu')
        ax.set_title("velocity map")
        ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
        ax.set_xticks(np.linspace(0, len(self.mmap), num=3))
        ax.set_xticklabels([str(-1*self.meta.rscale_lambda),"0", str(self.meta.rscale_lambda)])
        ax.set_ylabel("[kpc]")
        ax.set_yticks(np.linspace(0, len(self.mmap), num=5))
        yticks = ["{:.2f}".format(y) \
                    for y in np.linspace(-self.meta.reff*self.meta.rscale_lambda,
                                         self.meta.reff*self.meta.rscale_lambda, num=5)]
        ax.set_yticklabels(yticks)

#        ax.locator_params(tight=True, nbins=5)

    if do9plots:
# position and velocity histograms
        ax = axs[0,2]
        ax.hist(self.star['x'])
        ax.set_title("X pos")
        ax.locator_params(tight=True, nbins=5)

        ax = axs[1,2]
        ax.hist(self.star['y'])
        ax.set_title("Y pos")
        ax.locator_params(tight=True, nbins=5)

        ax = axs[2,0]
        ax.hist(self.star['vx'])
        ax.set_title("X vel")
        ax.locator_params(tight=True, nbins=5)

        ax = axs[2,1]
        ax.hist(self.star['vy'])
        ax.set_title("Y vel")
        ax.locator_params(tight=True, nbins=5)

    plt.tight_layout()
    if fn_save is None:
        fn_save = "galaxy_plot" + str(self.meta.id).zfill(5) + ".png"

    plt.savefig(fn_save, dpi=100)
    plt.close()

def save_gal(self, base='./'):

    def get_metadata(clazz):
        """
            Out of all attributes of a galaxy instance, leave only data.
        """
        return {name: attr for name, attr in clazz.__dict__.items()
                if not name.startswith("__")
                and not callable(attr)
                and not type(attr) is staticmethod}

    def get_metadata2(adict):
        return {name: attr for name, attr in adict.items()
                if not isinstance(attr, (np.ndarray, np.recarray, dict, list))}


    import h5py as hdf
    # Save data into a hdf5 file
    outfile = hdf.File(base + str(self.meta.id).zfill(6) + '_gal.hdf5',
                       'w', libver='latest')

    # Store metadata in HDF5 attributes
    attrs = get_metadata(self)
    attrs = get_metadata2(attrs)
    for name, atr in attrs.items():
        if atr != None:
            outfile.attrs.create(name, atr)
        #outfile.attrs[name] = atr

    # Store data under /selfaxy with direct assignment
    if hasattr(self, 'star'):
        #print("Saving star")
        star = outfile.create_group("star")
        for field in self.star.dtype.names:
            star.create_dataset(field, data=self.star[field])

    if hasattr(self, 'dm'):
        #print("Saving DM")
        dm = outfile.create_group("dm")
        for field in self.dm.dtype.names:
            dm.create_dataset(field, data=self.dm[field])

    if hasattr(self, 'cell'):
        #print("Saving gas")
        gas = outfile.create_group("gas")
        for field in self.cell.dtype.names:
            gas.create_dataset(field, data=self.cell[field])

    outfile.close()
