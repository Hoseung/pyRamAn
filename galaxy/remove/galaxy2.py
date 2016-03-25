# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 02:11:12 2015

@author: hoseung
"""

import numpy as np
#import  matplotlib.pyplot as plt
class Galaxy(object):

    def __init__(self, halo, radius_method='eff', info=None):
        self.halo = None
        self.set_halo(halo)
        self.id = int(halo['id'])
        self.radius_method=radius_method
        self.info = info
        self.xc = 0
        self.yc = 0
        self.zc = 0
        self.reff = 0
        self.mstar = 0
        self.nstar = 0
        self.mgas = 0
        self.lambda_arr = None
        self.lambda_r = 0
        self.vmap = None
        self.pt=[]
        self.pq=[]
        self.sfr=0.0
        self.ssfr=0.0

    def set_halo(self, halo):
        self.halo = halo
        
    def set_ptypes(self, pt=None):
        pass
   
    def physical_print(self, x):
        # If length, 
        
        print([i * self.info.pboxsize * 100 for i in list(x)])
    

    def mk_gal_from_gal(self, star, dm, cell,
               save=False, rscale=1.5, verbose=False,
               mstar_min=1e9,
               rmin = -1, Rgal_to_reff=5.0, method_com=2, follow_bp=None):
        """
            Input data in code unit.
            
            Parameters
            ----------
            
            Rgal_to_reff: 
                Galaxy radius = Reff * Rgal_to_reff.
                By default, Rgal = 5*Reff. 
                (E's have smaller R_t_r, L's have larger.)
                
            
        """
        import numpy as np
        member="Reff"
        print("Making a galaxy:", self.id)
        if verbose:
            print("SAVE:", save)
            print("RSCALE:", rscale)
            print("Halo size:", self.halo['rvir'] * self.info.pboxsize * 1000.0)

        rscale_cen = 0.25

        vx_hal0 = self.halo['vx']
        vy_hal0 = self.halo['vy']
        vz_hal0 = self.halo['vz']
        
        rr_tmp0 = min([self.halo['rvir'] * rscale_cen, 0.0002]) 
        # arbitrary! < 20kpc
        rr_tmp0 = max([rr_tmp0, 0.000025])
        # When merger occurs, larger radius is likely to include 
        # companion galaxy resulting center to be in the middle of nowhere.
        # If you want a larger galaxy, # increase rgal_tmp instead. 
        #        
        # xx is easier to search for than x.

        if verbose:
            print("First attempt to define the center"
                    " using particles inside the radius-0 from the halo center."
                    , rr_tmp0 * self.info.pboxsize * 1000, ['kpc'])
        
        xall = star['x']
        yall = star['y']
        zall = star['z']
        
        xc_tmp1 = np.median(xall)
        yc_tmp1 = np.median(yall)
        zc_tmp1 = np.median(zall)
        
        if verbose:
            print("Now guessed center is :", xc_tmp1, yc_tmp1, zc_tmp1)

        rr_tmp1 = min([self.halo['rvir'] * rscale, 0.00025]) # modified!
        if verbose:
            print("Searching for particles around the guessed center")
            print("within radius1", rr_tmp1 * self.info.pboxsize * 1000, ['kpc'])

        ind = np.where(np.square(np.subtract(xall, xc_tmp1)) \
                       + np.square(np.subtract(yall, yc_tmp1)) \
                       + np.square(np.subtract(zall, zc_tmp1)) < np.square(rr_tmp1))[0]

        xx = xall[ind]
        yy = yall[ind] 
        zz = zall[ind] 
        mm = star['m'][ind] 

        m_sum = sum(mm)
        print("Number of stars", len(ind))
        print("galaxy::mk_gal::m_sun  {:.3e}", m_sum * self.info.msun)
        if m_sum * self.info.msun < mstar_min:
            print("(2)Not enough stars: {:.2f} Msun".format(m_sum) * self.info.msun)
            print(" Aborting... \n")
            self.star = False
            return False

        # First, define the center w.r.t. stars.        

        if verbose:
            print("Now, cut out other components than the main galaxy \n" 
            " by identifying peaks in radial distribution of stars")
        if method_com == 1:
            if verbose: print("Second trial to determine the center. method1")
            xc0, yc0, zc0 = self.get_center(xx, yy, zz, mm) # STARS!
        elif method_com == 2:
            if verbose: print("Second trial to determine the center. method2")
            xc0, yc0, zc0 = self.get_center2(xx, yy, zz, mm) # STARS!
            
        
# Second, define the search radius
# Only 'simple' works, because 'star' is not defined yet.

        if verbose: 
            print("Currrent CENTER:", xc0,yc0,zc0)
            print("Calculate Effective radius")
        
        rr_tmp2 = self.reff_main_gal(xx-xc0, yy-yc0, zz-zc0) # No mass! 
        if rr_tmp2 == False:
            # no other galaxies
            rr_tmp2 = max([rr_tmp1, rmin])

        ind = np.where(np.square(xall - xc0) + np.square(yall - yc0) 
                            + np.square(zall - zc0) < np.square(rr_tmp2))[0]
        xx = np.subtract(xall[ind], xc0)
        yy = np.subtract(yall[ind], yc0)
        zz = np.subtract(zall[ind], zc0)
        mm = star['m'][ind]

        xc, yc, zc = self.get_center(xx, yy, zz, mm)
        xc = np.add(xc, xc0)
        yc = np.add(yc, yc0)
        zc = np.add(zc, zc0)
        
        if verbose: print("New CENTER:", xc,yc,zc)
        
        self.xc = xc
        self.yc = yc
        self.zc = zc

        ind = np.where(np.square(np.subtract(xall, xc)) \
                       + np.square(np.subtract(yall, yc)) \
                       + np.square(np.subtract(zall, zc)) < np.square(rr_tmp2))[0]
        
        if sum(star['m'][ind]) * self.info.msun < mstar_min:
            print("Not enough stars: {:.2f} Msun".format(sum(star['m'][ind]) * self.info.msun))
            print(" Aborting... \n")
            self.star = False
            return False                    

        xx = xall[ind] - xc
        yy = yall[ind] - yc
        zz = zall[ind] - zc
        mm = star['m'][ind]
        
        reff_tmp = self.get_radius(xx, yy, zz, mm, self.radius_method)
        rgal_tmp = min([rr_tmp2, reff_tmp * Rgal_to_reff])
        self.Rgal_to_reff = rgal_tmp / reff_tmp
        # should not surpass rr_tmp2, over where  another galaxy might be.
        
        # Test if the outder annulus has significant amount of stars
        # -> it shouldn't.        
        if verbose: print("Reff: ", reff_tmp * self.info.pboxsize*1000)

        if dm is not None:
            if member == "Reff":
                idm = np.where( np.square(dm["x"] - xc) + 
                                np.square(dm["y"] - yc) + 
                                np.square(dm["z"] - zc) <= np.square(rgal_tmp))[0]
            elif member == "v200":
            # Alghough the velocity is redefined later,
            # particle membership is fixed at this point. 
                idm = np.where( np.square(dm["vx"] - vx_hal0)+ 
                                np.square(dm["vy"] - vy_hal0)+ 
                                np.square(dm["vz"] - vz_hal0) <= np.square(200**2))[0]
            ndm_tot = len(idm)

            self.dm = np.recarray(ndm_tot, dtype=dm.dtype)
            if 'id' in dm.dtype.names:
                self.dm['id'] = dm['id'][idm] 
            if 'm' in dm.dtype.names:
                self.dm['m'] = dm['m'][idm] * self.info.msun
            if 'y' in dm.dtype.names:
                self.dm['x'] = (dm['x'][idm] - self.xc) * self.info.pboxsize * 1000
                self.dm['y'] = (dm['y'][idm] - self.yc) * self.info.pboxsize * 1000
                self.dm['z'] = (dm['z'][idm] - self.zc) * self.info.pboxsize * 1000
            if 'vy' in dm.dtype.names:
                # Velocity is centered and normalized later on.
                self.dm['vx'] = dm['vx'][idm]
                self.dm['vy'] = dm['vy'][idm]
                self.dm['vz'] = dm['vz'][idm]
            if verbose: print("DM data stored")
            
        if star is not None:
            istar = np.where(np.square(star["x"] - self.xc) + 
                            np.square(star["y"] - self.yc) + 
                            np.square(star["z"] - self.zc) <= np.square(rgal_tmp))[0]
            nstar_tot = len(istar)
            if verbose: print("nstar tot:", nstar_tot)        
            if verbose: print("Store stellar particle")
    
            self.star = np.recarray(nstar_tot, dtype=star.dtype)
            if 'id' in star.dtype.names:
                self.star['id'] = star['id'][istar]
            if 'm' in star.dtype.names:
                self.star['m'] = star['m'][istar] * self.info.msun
            if 'x' in star.dtype.names:
                self.star['x'] = (star['x'][istar] - self.xc) * self.info.pboxsize * 1000
            if 'y' in star.dtype.names:            
                self.star['y'] = (star['y'][istar] - self.yc) * self.info.pboxsize * 1000
            if 'z' in star.dtype.names:
                self.star['z'] = (star['z'][istar] - self.zc) * self.info.pboxsize * 1000
            if 'vx' in star.dtype.names:
                self.star['vx'] = star['vx'][istar]
            if 'vy' in star.dtype.names:                
                self.star['vy'] = star['vy'][istar]
            if 'vz' in star.dtype.names:                
                self.star['vz'] = star['vz'][istar]
            if 'time' in star.dtype.names:
                import utils.cosmology
                self.star['time'] = utils.cosmology.time2gyr(star['time'][istar],
                                             z_now = self.info.zred,
                                             info=self.info)
            if 'metal' in star.dtype.names:
                self.star['metal'] = star['metal'][istar]           

        if cell is not None:
            if verbose: print("Cell is NOT none")
            icell = np.where(np.square(cell["x"] - xc) + 
                np.square(cell["y"] - yc) + 
                np.square(cell["z"] - zc) <= np.square(rgal_tmp))[0]
            ncell_tot = len(icell)
            self.cell = np.recarray(ncell_tot, dtype=cell.dtype)
            self.cell['x'] = (cell['x'][icell] - self.xc) * self.info.pboxsize * 1000
            self.cell['y'] = (cell['y'][icell] - self.yc) * self.info.pboxsize * 1000
            self.cell['z'] = (cell['z'][icell] - self.zc) * self.info.pboxsize * 1000
            self.cell['dx'] = cell['dx'][icell]
            self.cell['var0'] = cell['var0'][icell]
            self.cell['var2'] = cell['var2'][icell]
            self.cell['var3'] = cell['var3'][icell]
            self.cell['var4'] = cell['var4'][icell]
            self.cal_mgas()
            
        # Some more sophistications.      
        
        #if self.reff < 0.5        
#        rgal_tmp = 4 * self.reff
        """
        print("Rgal = 4 * Reff = ", rgal_tmp * self.info.pboxsize * 1000)
        
            # Save sink particle as a BH, not cloud particles. 

        """
        self.reff = reff_tmp * self.info.pboxsize * 1000
        self.rgal = rgal_tmp * self.info.pboxsize * 1000
        self.nstar = nstar_tot
        self.mstar = sum(self.star['m'])
        #print(".........", self.star['m'][100:120], self.mstar)

        import utils.sampling as smp
        self.region = smp.set_region(xc=xc, yc=yc, zc=zc, radius = self.rgal)
        # Now, get cov        

        self.get_cov(center_only=True)
        self.star['vx'] = (self.star['vx'] - self.vxc) * self.info.kms
        self.star['vy'] = (self.star['vy'] - self.vyc) * self.info.kms
        self.star['vz'] = (self.star['vz'] - self.vzc) * self.info.kms
        self.dm['vx'] = (self.dm['vx'] - self.vxc) * self.info.kms
        self.dm['vy'] = (self.dm['vy'] - self.vyc) * self.info.kms
        self.dm['vz'] = (self.dm['vz'] - self.vzc) * self.info.kms

        return True

    def reff_main_gal(self, x, y, z):
        """
            Check for additional galaxy inside the region.
            If there is another galaxy, Rgal is limited to the 
            Finally returns temporary galaxy radius.
            
          1) peaks of components are far enough.
          2) minor component is larger than `10%' of the major component (number of stars)
           - second peak is also a considerable galaxy, large enough to be considered as
            a minor merging counter part. (or minor interaction counter part.)
          3) minimum value is significantly smaller than the second peak.
           - two galaxise are not closely connected.
           - minor components are not background from a larger structure. 
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
#        plt.ioff()
        #n, bins, patches = plt.hist(d, bins=nbins, normed=1)
        n, bins = np.histogram(d, bins=nbins, normed=1)

        imax = np.arange(len(n))[_find_maxima(_smooth(n,3), exclude_end=True)]
        # Because I ban the first bin being a maximum, the number of bins
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
            return self.halo['rvir']
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
        import numpy as np
        vdiff = np.sqrt(((vxt - vxt[iv])**2 +
                         (vyt - vyt[iv])**2 +
                         (vzt - vzt[iv])**2))
        return sum(-1./vdiff[vdiff != 0])
    

    def get_cov(self, multithreaded=False, center_only=False):
        import numpy as np

        if center_only:
            ind = np.where(self.star['x']**2 
                          +self.star['y']**2 
                          +self.star['z']**2 < self.reff )

            vx = self.star['vx'][ind]
            vy = self.star['vy'][ind]
            vz = self.star['vz'][ind]
            mm = self.star['m'][ind]
        else:
            vx = self.star['vx']
            vy = self.star['vy']
            vz = self.star['vz']

        vhalx = self.halo['vx']
        vhaly = self.halo['vy']
        vhalz = self.halo['vz']
        
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
            self.vxc = np.average(vxt[emin], weights=mm[emin]) # top 100 
            self.vyc = np.average(vyt[emin], weights=mm[emin])
            self.vzc = np.average(vzt[emin], weights=mm[emin])
        elif method == 'com':
            self.vxc, self.vyc, self.vzc = self.get_center(
                                vx, vy, vz, mm, tol=5, niter=8)

    def cal_mgas(self):
        msun = 1.98892e33 # solar mass in gram.
        self.mgas = sum((self.cell['var0'] * self.info.unit_d) * 
        (self.cell['dx'] * self.info.unit_l)**3 / msun)
        # [g/cm^3] * [cm^3] / [g/msun] = [msun]
    
    def cal_sfr(self):
        import utils
        import numpy as np
        # average over last 100Myrs
        time_new = 100 * 1e-3 # in Myr
        ind_new_star = np.where(utils.cosmology.time2myr(
                        self.star['time'], z_now=self.info.zred) < time_new)[0]
        m_star_new = sum(self.star['m'][ind_new_star])
        self.sfr=m_star_new / time_new
        return m_star_new / time_new
    
    def cal_ssfr(self):
        if self.sfr == 0:
            self.cal_sfr()
        self.ssfr = self.sfr / self.mstar

    def cal_vrot(self):
        import numpy as np
        d = np.sqrt(self.star['x']**2 +self.star['y']**2 +self.star['z']**2)
        dsort = np.sort(d)    
        pass

    def dist_map(self, npix=40):
        import numpy as np
        nx = npix
        ny = npix
        
        dist_map=np.zeros((nx, ny), dtype=float)
        for i in range(nx):
            for j in range(ny):
                dist_map[i][j]= np.sqrt((0.5 + i - nx/2)**2 + (0.5 + j - ny/2)**2)
        return dist_map

    def weighted_std(self, values, weights):
        import numpy as np
        import math
        """
        Return the weighted average and standard deviation.
    
        values, weights -- Numpy ndarrays with the same shape.
        """
        average = np.average(values, weights=weights)
        variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
           
        return math.sqrt(variance)
    
    def _photometry(self, img, nbin, npix_img, l_img):
        from Cappellari import mge
        import numpy as np
        
        arr=np.recarray((nbin,), dtype=[("eps", np.float),
                     ("mjr", np.float), ("pa", np.float),
                     ("xpos", np.float), ("ypos", np.float)])
        
    
        for i, frac in enumerate(np.linspace(0.005, 0.2, num=nbin, endpoint=True)):
            f = mge.find_galaxy.find_galaxy(img, quiet=True, plot=True, mask_shade=False, fraction=frac)
            arr.mjr[i] = f.majoraxis
            arr.eps[i] = f.eps
            arr.pa[i] = f.theta
            arr.xpos[i] = f.xpeak
            arr.ypos[i] = f.ypeak
    
        # convert units
        arr.mjr = arr.mjr * 3.5 / npix_img * l_img # in kpc.
        arr.xpos = (arr.xpos / npix_img - 0.5) * l_img # in kpc.
        arr.ypos = (arr.ypos / npix_img - 0.5) * l_img # in kpc.
            
        return arr
    
    def sigmap_vmap(self, npix, ind_ok=None):
        if ind_ok is None:
            rr = self.reff * self.rscale_lambda
    
            ind_ok = np.where((abs(self.star["x"]) <= rr) &
                              (abs(self.star["y"]) <= rr) &
                              (abs(self.star["z"]) <= rr) )[0]
            
            print(("\n" "Calculating Lambda_r for {} particles "
               "inside {:.3f}kpc, or {}Reff".format(len(ind_ok), rr, self.rscale_lambda)))
    
        x = self.star['x'][ind_ok]
        y = self.star['y'][ind_ok]
        m = self.star['m'][ind_ok]
        vz = self.star['vz'][ind_ok]
    
        # NGP charge assignment
        nx = npix
        ny = npix
    
        # distance in pixel unit.
        dist2d=np.zeros((nx, ny), dtype=float)
        for i in range(nx):
            for j in range(ny):
                dist2d[i][j]= np.sqrt((0.5 + i - nx/2)**2 + (0.5 + j - ny/2)**2)
    
        dist1d = dist2d.ravel()
    
        # fix center explicitely.
        # 0.5 * (min + max) != center
        xx = (x + rr) / rr * 0.5 * nx
        yy = (y + rr) / rr * 0.5 * ny
    
        ngx = np.clip(np.fix(xx), 0, nx-1)
        ngy = np.clip(np.fix(yy), 0, ny-1)
    
        indices = ngx + ngy * nx
    
        mmap = np.zeros(nx * ny, dtype=float)
        for i, ind in enumerate(indices):
            mmap[ind] += m[i]
    
        dx = (2 * rr) / npix
        mmap = mmap / (dx*dx)
        self.mmap = mmap.reshape(nx, ny)
    
        vmap = np.zeros(nx * ny, dtype=float)
        # mass-weighted sigma         
        sigmap=np.zeros(nx * ny, dtype=float)
        for i in range(nx * ny):
            ind = np.where(indices == i)[0]
            if len(ind) > 0:
                sigmap[i] = self.weighted_std(vz[ind], m[ind])
                vmap[i] = np.average(vz[ind], weights=m[ind])
            else:
                sigmap[i] = 0
                vmap[i] = 0
    
        self.sigmap = sigmap.reshape(nx, ny)
        self.vmap = vmap.reshape(nx, ny)
        #return mmap.reshape(nx, ny), sigmap.reshape(nx, ny), vmap.reshape(nx, ny)
        
        
    def general_ellipse(self, x, y, a, b, xoff, yoff, theta):
        return np.square((x - xoff) * np.cos(theta) + (y - yoff)*np.sin(theta)) / a**2 + \
        np.square((x - xoff) * np.sin(theta) - (y - yoff)*np.cos(theta)) / b**2
    
    def cal_lambda_better(self, npix_lambda=5,
                     r=0.5,
                     rscale=3.0,
                     method=1,
                     verbose=False):
        import numpy as np
        from Cappellari import mge
        import utils.sampling as smp
        from draw import pp

    
        # First, measure photometric properties
        nbin = int(npix_lambda * rscale) #  =15 
    
        npix_img = int(self.reff * rscale) * 4 # 5 * reff = Rgal in most case, 4 pixels in 1 kpc.
        self.npix_lambda = npix_img
        self.rscale_lambda = rscale
        region = smp.set_region(xc=0, yc=0, zc=0, radius = self.reff * rscale)  
        data = pp.den2d(self.star['x'], self.star['y'], self.star['z'], self.star['m'], \
              npix_img, region=region, cic=True, norm_integer=True)
    
        arr = self._photometry(data, nbin, npix_img, 2 * region['radius'])
        #arr.mjr = arr.mjr * 3.5 / npix_img *  # in kpc.
        
             
        # velocity map & velocity dispersion map.
        a_max = arr.mjr[-1] * 1.2 # 20% margin
        b_max = a_max * arr.eps[-1]
        xoff = np.mean(arr.xpos)
        yoff = np.mean(arr.ypos)
        #xoff = arr.xpos[-1]
        #yoff = arr.ypos[-1]
        theta = arr.pa[-1] # np.pi/2 - pa? 

               
        #ind_ok = np.where(q < 1)
#        print(npix_img)
        #mmap, sigmap, vmap = self.sigmap_vmap(npix_img, ind_ok=None)
        self.sigmap_vmap(npix_img, ind_ok=None)
     
#        x = x[ind_ok]
#        y = y[ind_ok]
#        m = self.star['m'][ind_ok]
#        vz = self.star['vz'][ind_ok]
        
        lambda_arr = np.zeros(nbin)
        # distance map        
        nx, ny = npix_img, npix_img
        if method == 1:
            xv, yv = np.meshgrid(np.arange(nx), np.arange(ny), sparse=False, indexing='ij')
            dist2d=np.zeros((nx, ny), dtype=float)
            for i in range(nx):
                for j in range(ny):
                    dist2d[i][j]= np.sqrt((0.5 + i - nx/2)**2 + (0.5 + j - ny/2)**2)
                    
       
        
        for i in range(nbin):
            a = arr.mjr[i]
            b = a * arr.eps[i]
            xoff = arr.xpos[i]
            yoff = arr.ypos[i]
            theta = np.pi/2 - arr.pa[i] #arr.pa[i] # np.pi/2 - pa? 

            if method == 1:
                q = self.general_ellipse(xv, yv, a, b, xoff, yoff, theta)
                ind = np.where((0.9 < q) & (q < 1.1))[0]  # 0.9 ? 1.1?
                a = np.sum(self.mmap[ind] * dist2d[ind] * abs(self.vmap[ind]), axis=None) # reduce 2D array into a scalar
                if a > 0:
                    ind2 = np.where(self.sigmap[ind] > 0)[0]
                    b = np.sum(self.mmap[ind[ind2]] * dist2d[ind[ind2]] 
                            * np.sqrt(self.vmap[ind[ind2]]**2 + self.sigmap[ind[ind2]]**2), axis=None)
                    lambda_arr[i] = a/b
            # particle
                
            elif method == 2:
                x = self.star['x'] - xoff                           
                y = self.star['y'] - yoff
                q = self.general_ellipse(x, y, a_max, b_max, 0, 0, theta)                        
                
                ind = np.where((0.9 < q) & (q < 1))[0] 
            
                if len(ind) > 0:
                    vv = abs(np.average(vz[ind], weights=m[ind]))
                    a = sum(m[ind] * dist[ind] * vv)
                    sig = self.weighted_std(vz[ind], m[ind])
                    b = sum(m[ind] * dist[ind] * np.sqrt(vv**2 + sig**2))
                    lambda_arr[i] = a/b
            
        self.lambda_arr = lambda_arr
        # npix = ind_1Reff. 
        # ** npix 2 = npix * rscale
        # 0.5 * npix = Reff
        # 0.25 * npix = 0.5Reff
        self.lambda_r = np.average(self.lambda_arr[int(0.2 * nbin) : int(0.2 * nbin) + 1])
        self.photometric = arr        
        if verbose: print("lambda_arr done")

        

  
        
    def cal_lambda_r(self, npix=10,
                     r=0.5,
                     rscale=3.0,
                     method=2,
                     verbose=False):
        """
        Parameters
        ----------
        npix: int
            number of points per 2*reff.
        r:
            ??
        rscale:
            Lambda_r is calculated for annuli up to Reff*rscale.
        
        Calculate rotation parameter(lambda_r).
        rotation parameter is calculated at all points (r < rscale*self.reff)
        constructing a curve. But only the value at 0.5Reff is used as 
        a representative value (Emsellem + 07).
        Values at r = 1.0 Reff is reported to behave similar to the value at 0.5Reff.
        
        method1:
        
        
        method2:
        
        """
        import numpy as np

        # already centered.
        self.rscale_lambda = rscale
        rr = self.reff * rscale
        npix2 = round(npix * rscale)

        ind_ok = np.where((abs(self.star["x"]) <= rr) &
                          (abs(self.star["y"]) <= rr) &
                          (abs(self.star["z"]) <= rr) )[0]
        x = self.star['x'][ind_ok]
        y = self.star['y'][ind_ok]
        m = self.star['m'][ind_ok]
        vz = self.star['vz'][ind_ok]
# Some margin makes mmap and vmap look better.
# But only inside Lambda_R for R < Reff is returned. 
        print(("\n" "Calculating Lambda_r for {} particles "
               "inside {:.3f}kpc, or {}Reff".format(len(ind_ok), rr, rscale)))
        
        # NGP charge assignment
        nx = npix2
        ny = npix2

        # distance in pixel unit.
        dist2d=np.zeros((nx, ny), dtype=float)
        for i in range(nx):
            for j in range(ny):
                dist2d[i][j]= np.sqrt((0.5 + i - nx/2)**2 + (0.5 + j - ny/2)**2)
                
        dist1d = dist2d.ravel()

        # fix center explicitely.
        # 0.5 * (min + max) != center
        xx = (x + rr) / rr * 0.5 * nx
        yy = (y + rr) / rr * 0.5 * ny
        
        ngx = np.clip(np.fix(xx), 0, nx-1)
        ngy = np.clip(np.fix(yy), 0, ny-1)
        
        indices = ngx + ngy * nx

        mmap = np.zeros(nx * ny, dtype=float)
        for i, ind in enumerate(indices):
            mmap[ind] += m[i]
        
        dx = (2 * rr) / npix2
        mmap = mmap / (dx*dx)
        self.mmap = mmap.reshape(nx, ny)
        
        vmap = np.zeros(nx * ny, dtype=float)
        # mass-weighted sigma         
        sigmap=np.zeros(nx * ny, dtype=float)
        for i in range(nx * ny):
            ind = np.where(indices == i)[0]
            if len(ind) > 0:
                sigmap[i] = self.weighted_std(vz[ind], m[ind])
                vmap[i] = np.average(vz[ind], weights=m[ind])
            else:
                sigmap[i] = 0
                vmap[i] = 0

        self.dist1d = dist1d
        self.sigmap = sigmap.reshape(nx, ny)
        self.vmap = vmap.reshape(nx,ny)

        points = np.zeros(0.5 * npix2) # 1.5 Reff
        npoints = len(points)
        if method == 1:
            for i in range(npoints):
                ind = np.where((dist1d > i) & (dist1d < i+1))[0]
                a = sum(mmap[ind] * dist1d[ind] * abs(vmap[ind]))
                if a != 0:
                    ind2 = np.where(sigmap[ind] > 0)[0]
                    b = sum(mmap[ind[ind2]] * dist1d[ind[ind2]] 
                            * np.sqrt(vmap[ind[ind2]]**2 + sigmap[ind[ind2]]**2))
                    points[i] = a/b
        elif method == 2:
            dd = np.sqrt(x**2 + y**2)
            for i in range(npoints):
                ind = np.where( (dd > i*rr/npoints) & (dd < (i+1)*rr/npoints))[0]
                if len(ind) > 0:
                    vv = abs(np.average(vz[ind], weights=m[ind]))
                    a = sum(m[ind] * dd[ind] * vv)
                    sig = self.weighted_std(vz[ind], m[ind])
                    b = sum(m[ind] * dd[ind] * np.sqrt(vv**2 + sig**2))
                    points[i] = a/b
            
        self.lambda_arr = points
        # npix = ind_1Reff. 
        # ** npix 2 = npix * rscale
        # 0.5 * npix = Reff
        # 0.25 * npix = 0.5Reff
        self.lambda_r = np.average(self.lambda_arr[int(0.25 * npix) : int(0.25 * npix) + 1])
        if verbose: print("lambda_arr done")
#                                + self.lambda_arr[len(self.lambda_arr) -1])
    
    def reorient(self, part, modify=False):
        """
        part can be of any type.
        It's stellar particles by default.
        """
#        import numpy as np
        if not hasattr(self, 'nvec'):
            self.cal_norm_vec(part)

        try:
            RM = self.rotation_matrix
        except:
            print("Galaxy.reorientation: No Rotation matrix. \n Calculating..")
            RM = self.cal_rotation_matrix()       
        
        # New coordinates
        # how about matrix dot operation?
        new_x = RM[0,0]*part['x'] + RM[1,0]*part['y'] + RM[2,0]*part['z']
        new_y = RM[0,1]*part['x'] + RM[1,1]*part['y'] + RM[2,1]*part['z']
        new_z = RM[0,2]*part['x'] + RM[1,2]*part['y'] + RM[2,2]*part['z']

        new_vx= RM[0,0]*part['vx'] + RM[1,0]*part['vy'] + RM[2,0]*part['vz']
        new_vy= RM[0,1]*part['vx'] + RM[1,1]*part['vy'] + RM[2,1]*part['vz']
        new_vz= RM[0,2]*part['vx'] + RM[1,2]*part['vy'] + RM[2,2]*part['vz']
        
        if modify:
            part['x'] = new_x
            part['y'] = new_y
            part['z'] = new_z
            part['vx'] = new_vx
            part['vy'] = new_vy
            part['vz'] = new_vz
        else:
            return new_x, new_y, new_z

    def cal_norm_vec(self, part=None):
        """
        Determines normal vector of the given galaxy.
        'galaxy' means stellar particles by default. 
        Add self.nvec attribute and also returns nvec
        """
        import numpy as np
        if part is None:
            part = self.star

        x = part['x']
        y = part['y']
        z = part['z']
        vx= part['vx']
        vy= part['vy']
        vz= part['vz']

        lx= sum(y*vz - z*vy) # normalized; *m is omitted.
        ly= sum(z*vx - x*vz)
        lz= sum(x*vy - y*vx)

        self.nvec = np.array([lx, ly, lz])/np.sqrt(lx**2 + ly**2 + lz**2)
        return self.nvec


    def _get_rotation_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        http://stackoverflow.com/a/6802723
        """
        import numpy as np
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
        import numpy as np
        import numpy.linalg as lag
        import math
        # rotation axis
        dest = np.asarray(dest) # towards +z direction.
        if nvec is None:
            try:
                nvec = self.nvec
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
            Follow 
        
        """
        
        
    def bound_particles(self, n=100, ptype='star'):
        """
            save IDs of n most bound particles as galaxy.most_bound_particles.
            
            NOTE
            ----
                Requires center of mass and center of velocity already determined.
        """
        
        #vrange = 100
        vx = self.star['vx']
        vy = self.star['vy']
        vz = self.star['vz']
    
        vdiff = np.square(vx - self.vxc) + np.square(vy - self.vyc) + np.square(vz - self.vzc)
        npart = len(vdiff)
        
        if npart < 100:
            print("Too little particles within velocity window km/s")
            return False
    
        vdiff = np.zeros(npart, dtype=float)
        
        emin = np.argsort(vdiff)[:min(n, npart)]
        self.most_bound_particles = self.star['id'][emin]
        
    
    def cal_b2t(self, ptype='star', hist=False):
        """
        Reorient galaxy first. 
        
        parameters
        ----------------------
        hist : 
            if True, ellipticity histogram is returned
        """
        import numpy as np
        # Radius

        # Calculating boundness requires total mass inside a radius.
        # -> DM, Cell are also needed. 
        part = getattr(self, ptype)
        x = part['x']
        y = self.star['y']
        z = self.star['z']

        m_all_arr = np.concatenate((self.star['m'], 
                              self.dm['m'], 
                              self.cell['var0']*self.info.unit_d
                              *(self.cell['dx']*self.info.unit_l)**3/1.98e33))

        rstar = np.sqrt(np.array(x**2 + y**2 + z**2))

        r_all_arr = np.concatenate((rstar,
                                np.sqrt(self.dm['x']**2 
                                + self.dm['y']**2 
                                + self.dm['z']**2),
                                np.sqrt(self.cell['x']**2
                                + self.cell['y']**2
                                + self.cell['z']**2)))

        ind_r_all_sort = np.argsort(r_all_arr)
        m_all_radial = np.cumsum(m_all_arr[ind_r_all_sort])
        
        # enclosed mass at the location of each star particle.
        enclosed_mass_at_r = np.interp(rstar, r_all_arr[ind_r_all_sort], m_all_radial)
        
        boxtokpc = self.info.pboxsize * 1000
        G = 6.67384e-11  # m^3 kg^-1 s^-2
        kpc_to_m = 3.08567758e19 
        msun = 1.9891e30 # kg
        v_circ = np.sqrt(G * msun * enclosed_mass_at_r 
                        /(kpc_to_m * rstar * boxtokpc)) * 1e-3 # m/s to in km/s

        j_circ = rstar * v_circ # * boxtokpc omitted (later compensated.)
        # Finally, r in code unit, v in km/s

        j_phi = (x * self.star['vy'] - y * self.star['vx']) # * boxtokpc omitted.
        ellipticity = j_phi / j_circ
                   
        # < 0.8 from Scannapieco 2009
        disk = sum(self.star['m'][np.where(ellipticity < 0.8)])
        self.d2t = disk / self.mstar # or sum(self.star['m'][ind])
        if hist:
            return np.histogram(ellipticity, range=[-2,2], bins=20)

#    def _radial_profile(x,y,m):
#        dsored = np.sqrt(x**2 + y**2)
#        inds = np.d
     

    def cal_radial_profile(self):
        x = self.star['x']
        y = self.star['y']
        m = self.star['m']
        
        # x-y plane
        f = _radial_mass(x,y,m)           
        

    def cal_trivia(self):
        self.mstar = sum(self.star['m'])
        self.vcen = self.get_vcen()
#        self.metal = 
#        self.reff = 12333

    def save_gal_pickle(self):
        import pickle
        with open(str(self.id).zfill(6) + 'gal.pickle', 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
         
    def plot_gal(self, fn_save=None, ioff=True, **kwargs):
        import matplotlib.pyplot as plt
        from matplotlib import ticker
        from draw import pp
        from matplotlib.colors import LogNorm
#        import utils.prettyplot as ptt
        import numpy as np
        if ioff:
            plt.ioff()
            
        npix = self.npix_lambda

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


        rgal_to_reff = self.Rgal_to_reff

        do9plots = False        
        if do9plots:
            fig, axs = plt.subplots(3,3)
            fig.set_size_inches(12,9) # 9 plots
        else:
            fig, axs = plt.subplots(2,2)

        fig.suptitle("ID: {}    z: {:2f}".format(str(self.id).zfill(5), self.info.zred))

# Stellar particle density map
        ax = axs[0,0]
        """
            1st plot. - Stellar density map. 
        """        
        img = pp.part2den(self.star, self.info, npix=npix)
#        region=self.region, offset=[self.xc, self.yc, self.zc])

        im = ax.imshow(img.data, cmap=plt.get_cmap('brg'),
                       norm=LogNorm(vmin=1e6))
        cb = plt.colorbar(im, ax=ax)

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
                    for y in np.linspace(-self.rgal, self.rgal, num=5)]
        ax.set_yticklabels(yticks)

# Lambda_r plot
        if self.lambda_arr is not None:
            ax = axs[0,1]
            ax.plot(self.lambda_arr) # ~ 1 * Reff
            ax.set_title(r"$\lambda _{R}$") 
            ax.text(0.5 * len(self.lambda_arr), 0.8, "{:.2e}".format(self.mstar) + r"$M_{\odot}$")
            ax.set_ylim(bottom=0, top=1)
            ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
            
            ll = self.rscale_lambda
            unit = len(self.lambda_arr)/ll # number of points per 1Reff
            xticks = [0, 0.5*unit]
            reffs=[0, 0.5]
            for i in range(int(ll)): # If lambda_r is calculated furthrer -
                xticks.append(unit * (i+1))
                reffs.append((i+1))
            ax.set_xticks(xticks)
            xticks_label = ["{:.1f}".format(rf) for rf in reffs]
            ax.set_xticklabels(xticks_label)

# sigma map
        if hasattr(self, "sigmap"):
            l_range = self.reff * self.rscale_lambda
            ax = axs[1,0]        
            im = ax.imshow(self.sigmap, vmin=0, vmax = 200, cmap=plt.get_cmap('brg'))
            cb = plt.colorbar(im, ax=ax, label=r'$km s^{-1}$') # needs a mappable object.
            tick_locator = ticker.MaxNLocator(nbins=3)
            cb.locator = tick_locator
            cb.update_ticks()
            ax.set_title(r"$\sigma$ map")
            ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
            ax.set_xticks(np.linspace(0, len(self.mmap), num=3))
            ax.set_xticklabels([str(-1*self.rscale_lambda),"0", str(self.rscale_lambda)])
            ax.set_ylabel("[kpc]")
            ax.set_yticks(np.linspace(0, len(self.mmap), num=5))
            yticks = ["{:.2f}".format(y) \
                        for y in np.linspace(-l_range, l_range, num=5)]
            ax.set_yticklabels(yticks)        
#        ax.locator_params(tight=True, nbins=5)

# velocity map
        if self.vmap is not None:
            l_range = self.reff * self.rscale_lambda
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
            ax.set_title("velocity map")
            ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
            ax.set_xticks(np.linspace(0, len(self.mmap), num=3))
            ax.set_xticklabels([str(-1*self.rscale_lambda),"0", str(self.rscale_lambda)])
            ax.set_ylabel("[kpc]")
            ax.set_yticks(np.linspace(0, len(self.mmap), num=5))
            yticks = ["{:.2f}".format(y) \
                        for y in np.linspace(-self.reff*self.rscale_lambda,
                                             self.reff*self.rscale_lambda, num=5)]
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
            fn_save = "galaxy_plot" + str(self.id).zfill(5) + ".png"

#        from os.path import isdir
#        from os import mkdir
#        if not isdir(save_dir):
#            mkdir(save_dir)
        plt.savefig(fn_save, dpi=200)
        plt.close()
