# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 23:45:55 2015

Measure ICL fraction by merger channel

@author: hoseung
"""
# Galaxy object 
#class Galaxy(halo.HaloMeta):
"""
Inherites halo.HaloMeta. HaloMeta does not contain 
positional and kinematic information of the halo. This is good because
recalculated center of galaxy is likely to deviate from the center of the halo. 
"""
#class Galaxymeta():
    
import numpy as np
#import  matplotlib.pyplot as plt
class Galaxy(object):

    def __init__(self, halo=None, radius_method='eff', info=None):
        self.set_halo(halo)
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
        self.mrj=0.0
        self.d2t=0.0

    def set_halo(self, halo):
        if halo is not None:
            self.halo = halo
            self.id = int(halo['id'])
        
    def set_ptypes(self, pt=None):
        pass
   
    def physical_print(self, x):
        # If length, 
        
        print([i * self.info.pboxsize * 100 for i in list(x)])
    
    def mk_gal(self, star, dm, cell,
               save=False, rscale=0.5, verbose=False,
               mstar_min=1e9,
               rmin = -1, Rgal_to_reff=5.0, method_com=1, follow_bp=None):
        """
            Input data in code unit.
            
            Parameters
            ----------
            
            Rgal_to_reff: 
                Galaxy radius = Reff * Rgal_to_reff.
                By default, Rgal = 5*Reff. 
                (E's have smaller R_t_r, L's have larger.)
                
        """
        member="Reff"
        print("Making a galaxy:", self.id)
        if verbose:
            print("SAVE:", save)
            print("RSCALE:", rscale)
            print("Halo size:", self.halo['rvir'] * self.info.pboxsize * 1000.0)

        rscale_cen = 0.25
        #xc_tmp0 = self.halo['x']
        #yc_tmp0 = self.halo['y']
        #zc_tmp0 = self.halo['z']
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

        xall = star['pos'][:,0]
        yall = star['pos'][:,1]
        zall = star['pos'][:,2]
        
        xc_tmp1 = np.median(xall)
        yc_tmp1 = np.median(yall)
        zc_tmp1 = np.median(zall)
        if verbose:
            print("Now guessed center is :", xc_tmp1, yc_tmp1, zc_tmp1)

        rr_tmp1 = min([self.halo['rvir'] * rscale, 0.00025]) # modified!
        if verbose:
            print("Search for particles around the guessed center")
            print("within radius1", rr_tmp1 * self.info.pboxsize * 1000, ['kpc'])

        ind = np.where(np.square(xall - xc_tmp1) + np.square(yall - yc_tmp1) 
                            + np.square(zall - zc_tmp1) < np.square(rr_tmp1))[0]

        xx = xall[ind]
        yy = yall[ind] 
        zz = zall[ind] 
        mm = star['m'][ind] 

        m_sum = sum(mm)
        print("Number of stars", len(ind))
        print("galaxy::mk_gal: {:.2e} ".format(m_sum * self.info.msun))
        if m_sum * self.info.msun < mstar_min:
            print("(2)Not enough stars")
            print(" Aborting... \n")
            self.star = False
            return False
#        else:
            #if verbose:
#            print("\n Galaxy is made. \n Number of particles:", len(ind))

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
        #rr_tmp2 = self.reff_main_gal(xall-xc0, yall-yc0, zall-zc0) # No mass! 
        rr_tmp2 = self.reff_main_gal(xx-xc0, yy-yc0, zz-zc0) # No mass! 
        if rr_tmp2 == False:
            # no other galaxies
            rr_tmp2 = max([rr_tmp1, rmin])

        ind = np.where(np.square(xall - xc0) + np.square(yall - yc0) 
                            + np.square(zall - zc0) < np.square(rr_tmp2))[0]
        xx = xall[ind] - xc0
        yy = yall[ind] - yc0
        zz = zall[ind] - zc0
        mm = star['m'][ind]

        xc, yc, zc = self.get_center(xx, yy, zz, mm)
        xc = xc + xc0
        yc = yc + yc0
        zc = zc + zc0
        if verbose: print("New CENTER:", xc,yc,zc)
        self.xc = xc
        self.yc = yc
        self.zc = zc

        ind = np.where(np.square(xall - xc) + np.square(yall - yc)
                     + np.square(zall - zc) < np.square(rr_tmp2))[0]
        
        if sum(star['m'][ind]) * self.info.msun < mstar_min:
            print("Not enough stars: {:.2e} Msun".format(sum(star['m'][ind]) * self.info.msun))
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
                idm = np.where( np.square(dm["pos"][:,0] - xc) + 
                                np.square(dm["pos"][:,1] - yc) + 
                                np.square(dm["pos"][:,2] - zc) <= np.square(rgal_tmp))[0]
            elif member == "v200":
            # Alghough the velocity is redefined later,
            # particle membership is fixed at this point. 
                idm = np.where( np.square(dm["pos"][:,0] - vx_hal0)+ 
                                np.square(dm["pos"][:,1] - vy_hal0)+ 
                                np.square(dm["pos"][:,2] - vz_hal0) <= np.square(200**2))[0]
            ndm_tot = len(idm)

            self.dm = np.recarray(ndm_tot, dtype=dm.dtype)
            if 'id' in dm.dtype.names:
                self.dm['id'] = dm['id'][idm] 
            if 'm' in dm.dtype.names:
                self.dm['m'] = dm['m'][idm] * self.info.msun
            if 'pos' in dm.dtype.names:
                self.dm['pos'] = (dm['pos'][idm,:] - [self.xc, self.yc, self.zc]) * self.info.pboxsize * 1000
#                self.dm['pos'][:,1] = (dm['pos'][idm,1] - self.yc) * self.info.pboxsize * 1000
#                self.dm['pos'][:,2] = (dm['pos'][idm,2] - self.zc) * self.info.pboxsize * 1000
#                self.dm['y'] = (dm['y'][idm] - self.yc) * self.info.pboxsize * 1000
#                self.dm['z'] = (dm['z'][idm] - self.zc) * self.info.pboxsize * 1000
            if 'vel' in dm.dtype.names:
                # Velocity is centered and normalized later on.
                self.dm['vel'] = dm['vel'][idm,:]
#                self.dm['vy'] = dm['vy'][idm]
#                self.dm['vz'] = dm['vz'][idm]
            if verbose: print("DM data stored")
            
        if star is not None:
            istar = np.where(np.square(star["pos"][:,0] - self.xc) + 
                             np.square(star["pos"][:,1] - self.yc) + 
                             np.square(star["pos"][:,2] - self.zc) <= np.square(rgal_tmp))[0]
            nstar_tot = len(istar)
            if verbose: print("nstar tot:", nstar_tot)        
            if verbose: print("Store stellar particle")
    
            self.star = np.recarray(nstar_tot, dtype=star.dtype)
            if 'id' in star.dtype.names:
                self.star['id'] = star['id'][istar]
            if 'm' in star.dtype.names:
                self.star['m'] = star['m'][istar] * self.info.msun
            if 'pos' in star.dtype.names:
                self.star['pos'] = (star['pos'][istar,:] - [self.xc, self.yc, self.zc]) * self.info.pboxsize * 1000
#                self.star['pos'][:,1] = (star['pos'][istar,1] - self.xc) * self.info.pboxsize * 1000
#                self.star['pos'][:,2] = (star['pos'][istar,2] - self.xc) * self.info.pboxsize * 1000
#            if 'y' in star.dtype.names:            
#                self.star['y'] = (star['y'][istar] - self.yc) * self.info.pboxsize * 1000
#            if 'z' in star.dtype.names:
#                self.star['z'] = (star['z'][istar] - self.zc) * self.info.pboxsize * 1000
            if 'vel' in star.dtype.names:
                self.star['vel'] = star['vel'][istar,:]
#            if 'vy' in star.dtype.names:                
#                self.star['vy'] = star['vy'][istar]
#            if 'vz' in star.dtype.names:                
#                self.star['vz'] = star['vz'][istar]
            if 'time' in star.dtype.names:
                import utils.cosmology
                self.star['time'] = utils.cosmology.time2gyr(star['time'][istar],
                                             z_now = self.info.zred,
                                             info=self.info)
            if 'metal' in star.dtype.names:
                self.star['metal'] = star['metal'][istar]           

        if cell is not None:
            if verbose: print("Cell is NOT none")
            icell = np.where(np.square(cell["pos"][:,0] - xc) + 
                             np.square(cell["pos"][:,1] - yc) + 
                             np.square(cell["pos"][:,2] - zc) <= np.square(rgal_tmp))[0]
            ncell_tot = len(icell)
            dtype_cell = [('pos', '<f8', (3,)), ('dx', '<f8'),
                          ('rho', '<f8'), ('vel', '<f8', (3,)),
                          ('temp', '<f8'), ('metal', '<f8')]

            self.cell = np.recarray(ncell_tot, dtype=dtype_cell)
            self.cell['pos'] = (cell['pos'][icell,:] - [self.xc, self.yc, self.zc]) * self.info.pboxsize * 1000
#            self.cell['pos'][:,1] = (cell['pos'][icell,1] - self.xc) * self.info.pboxsize * 1000
#            self.cell['pos'][:,2] = (cell['pos'][icell,2] - self.xc) * self.info.pboxsize * 1000

            #self.cell['y'] = (cell['y'][icell] - self.yc) * self.info.pboxsize * 1000
            #self.cell['z'] = (cell['z'][icell] - self.zc) * self.info.pboxsize * 1000
            self.cell['dx'] = cell['dx'][icell]
            self.cell['rho'] = cell['var0'][icell]
            self.cell['vel'][:,0] = cell['var1'][icell]
            self.cell['vel'][:,1] = cell['var2'][icell]
            self.cell['vel'][:,2] = cell['var3'][icell]
            self.cell['temp'] = cell['var4'][icell]
            self.cell['metal'] = cell['var5'][icell]
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
        self.star['vel'] = (self.star['vel'] -[self.vxc, self.vyc, self.vzc]) * self.info.kms
        #self.star['vy'] = (self.star['vy'] - self.vyc) * self.info.kms
        #self.star['vz'] = (self.star['vz'] - self.vzc) * self.info.kms
        self.dm['vel'] = (self.dm['vel'] - [self.vxc,self.vyc,self.vzc]) * self.info.kms
        #self.dm['vy'] = (self.dm['vy'] - self.vyc) * self.info.kms
        #self.dm['vz'] = (self.dm['vz'] - self.vzc) * self.info.kms

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
        self.mgas = sum((self.cell['rho'] * self.info.unit_d) * 
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


    def _pseudo_particles(self,xx,yy,mm,vz, sig=0.5, n_times=60):
            cov = [[sig,0],[0,sig]]
            dx, dy = np.random.multivariate_normal([0,0], cov, len(xx) * n_times).T
            xnew = np.repeat(xx, n_times) + dx
            ynew = np.repeat(yy, n_times) + dy
            mm = np.repeat(mm,n_times)/n_times
            vz = np.repeat(vz,n_times)
            return xnew, ynew, mm, vz

    def cal_lambda_r_eps(self,
                         npix_per_reff=5,
                         rscale=None,
                         method='ellip',
                         verbose=False,
                         galaxy_plot_dir='./',
                         n_pseudo=1,
                         voronoi=None):
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
        
        Fixed number of pixels for all galaxies.
        boxsize = rscale * Reff * 2
        dx = Reff / npix
        npix_all = npix * rscale
        
        method1:
        
        
        method2:
        
        """
        import numpy as np

        # already centered.
        if rscale is not None:
            self.rscale_lambda = rscale
            print("self.rscale_lambda is updated to {}.".format(rscale))
        else:
            rscale = self.rscale_lambda
            
        if rscale is None:
            print("Error!")
            return


        if n_pseudo > 1:
            # sig in kpc unit.
            n_pseudo = max([round(1e6/self.nstar), n_pseudo])
            xstars, ystars, mm, vz = self._pseudo_particles(self.star['x'],
                                                      self.star['y'],
                                                      self.star['m'],
                                                      self.star['vz'],
                                                      sig=0.3,
                                                      n_times = n_pseudo)
        else:
            xstars = self.star['x']
            ystars = self.star['y']
            mm = self.star['m']
            vz = self.star['vz']

        reff = self.reff    
        rr = reff * rscale
        dx = reff / npix_per_reff
        
        # Some margin makes mmap and vmap look better.
        # If rscale = 3 is given, calculate inside 3Reff,
        # but plot a map of 4Reff.
        l_img = 2 * (rscale + 1) * reff
        npix = round(npix_per_reff * 2 * (rscale + 1)) 
        points = np.zeros(round(npix_per_reff * rscale))
        npoints = len(points)

        nx, ny = npix, npix
    
        if verbose: print(("\n" "Calculating Lambda_r for {} particles " 
        "inside {:.3f}kpc, or {}Reff".format(len(self.star['x']), rr, rscale)))
        
        
    # Construct mass map, velocity map, sigma map.
        # using NGP charge assignment

        # fix center explicitely.
        # 0.5 * (min + max) != center

        xstars = (xstars + rr) / rr * 0.5 * nx
        ystars = (ystars + rr) / rr * 0.5 * ny
        # because of Gaussian smoothing, pseudo particles can go out of cic region.
        # But don't worry, np.clip is ready.
        
        ngx = np.clip(np.fix(xstars), 0, nx-1) 
        ngy = np.clip(np.fix(ystars), 0, ny-1)
        
        indices = (ngx + ngy * nx).astype(np.int32)

        mmap = np.zeros(nx * ny, dtype=float)

        if voronoi is not None:
            count_map = np.zeros(nx * ny, dtype=float)
            for i, ind in enumerate(indices):
                count_map[ind] += 1
                mmap[ind] += mm[i]
        else:
            for i, ind in enumerate(indices):
                mmap[ind] += mm[i]

        mmap = mmap / (dx*dx)
        self.mmap = mmap.reshape(nx, ny)
        
        vmap = np.zeros(nx * ny, dtype=float)
        # mass-weighted sigma         
        sigmap=np.zeros(nx * ny, dtype=float)


        if voronoi is not None:
            i_ok = np.unique(indices) # non-zero only?
#            print(i_ok)
            count_map = count_map[i_ok]
            noise_map = np.sqrt(count_map)
 
            xpos_regular = np.repeat(np.arange(nx),ny)[i_ok]
            ypos_regular = np.tile(np.arange(ny),nx)[i_ok]
            from Cappellari.voronoi.voronoi_2d_binning import voronoi_2d_binning
            """
            This function accepts only data on uniform grid...?
            """
            binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = \
                voronoi_2d_binning(xpos_regular, ypos_regular, count_map,
                                 noise_map, targetSN=voronoi["targetSN"],
                                 plot=voronoi["plot"], quiet=voronoi["quiet"])



            def _display_pixels(x, y, counts, pixelSize):
                """
                Display pixels at coordinates (x, y) coloured with "counts".
                This routine is fast but not fully general as it assumes the spaxels
                are on a regular grid. This needs not be the case for Voronoi binning.
            
                """
                xmin, xmax = np.min(x), np.max(x)
                ymin, ymax = np.min(y), np.max(y)
                nx = round((xmax - xmin)/pixelSize) + 1
                ny = round((ymax - ymin)/pixelSize) + 1
                img = np.full((nx, ny), np.nan)  # use nan for missing data
                j = np.round((x - xmin)/pixelSize).astype(int)
                k = np.round((y - ymin)/pixelSize).astype(int)
                img[j, k] = counts
            
                plt.imshow(np.rot90(img), interpolation='none', cmap='prism',
                           extent=[xmin - pixelSize/2, xmax + pixelSize/2,
                                   ymin - pixelSize/2, ymax + pixelSize/2])
           
           
                                                         
            mmap_v = np.zeros(len(xNode))
            vmap_v = np.zeros(len(xNode)) # Not actually maps, but 1-D arrays.
            sigmap_v = np.zeros(len(xNode))

            # Quantities in Voronoi bin
            for ibin in np.arange(len(xNode)):
                ind = np.where(binNum == ibin)[0] # pixels in this Voronoi bin
                i_part = np.empty((0,0), dtype=int)
                for j in ind:
                    i_part = np.append(i_part, np.where(indices == j)[0])
                # all particles belonging to one Voronoi cell
                mmap_v[ibin] = sum(mm[i_part])
                try:
                    vmap_v[ibin] = np.average(vz[i_part], weights=mm[i_part])
                except:
                     continue # 
 #                   print("mass", mm[i_part])
#                    break
                sigmap_v[ibin] = self.weighted_std(vz[i_part], weights=mm[i_part])

                # update original map too.
                vmap[ind] = vmap_v[ibin]
                sigmap[ind] = sigmap_v[ibin]

            lambdamap_v = vmap_v / sigmap_v


            mmap_org = mmap 
            vmap_org = vmap 
            sigmap_org = sigmap 
            self.sigmap = sigmap_org.reshape(nx, ny)
            self.vmap = vmap_org.reshape(nx,ny)

            # npix * npix map for plots
            for ibin in range(len(xNode)):
                ind = np.where(binNum == ibin)[0]
                vmap_org[ind] = vmap_v[ibin]
                sigmap_org[ind] = sigmap_v[ibin]

            import matplotlib.pyplot as plt
            
            rnd = np.argsort(np.random.random(xNode.size))  # Randomize bin colors
            _display_pixels(xpos_regular, ypos_regular, rnd[binNum], nx)
#            plt.imshow(self.vmap)
            plt.savefig('./mmap.png')
            plt.imshow(self.sigmap)
            plt.savefig('./sig.png')
            plt.close()

            # update maps.
            # keep oroginal mmap to draw plots
            # Don't worry.  it's just switching pointer. No memory copy.
            sigmap = sigmap_v
            vmap = vmap_v # overwrite??
            mmap = mmap_v          

        else:
            for i in range(nx * ny):
                ind = np.where(indices == i)[0]
                if len(ind) > 0:
                    sigmap[i] = self.weighted_std(vz[ind], mm[ind])
                    vmap[i] = np.average(vz[ind], weights=mm[ind])
                else:
                    sigmap[i] = 0
                    vmap[i] = 0
            xNode = np.tile(np.arange(nx),ny) # x = tile? or repeat? 
            yNode = np.repeat(np.arange(ny),nx)
            self.sigmap = sigmap.reshape(nx, ny)
            self.vmap = vmap.reshape(nx,ny)

# calculate profile by summing in bins.
        if method == 'circ':
            xcen = nx/2
            ycen = ny/2
            dist1d = np.sqrt(np.square(xNode - xcen) + np.square(yNode - ycen))

            for i in range(npoints):
                ind = np.where((dist1d > i) & (dist1d < (i+1)))[0]
#                print(i*rr/npoints, len(ind))
                a = sum(mmap[ind] * dist1d[ind] * abs(vmap[ind]))
                if a != 0:
                    ind2 = np.where(sigmap[ind] > 0)[0]
                    b = sum(mmap[ind[ind2]] * dist1d[ind[ind2]] 
                            * np.sqrt(vmap[ind[ind2]]**2 + sigmap[ind[ind2]]**2))
                    points[i] = a/b

            if voronoi:
                dd = np.sqrt( (xNode - 0.5*npix)**2 + (yNode - 0.5*npix)**2 )
                i_radius = np.fix(dd).astype(int)
                new_arr = np.zeros(np.fix(max(dd)) + 1)
                new_cnt = np.zeros(np.fix(max(dd)) + 1)
                # NGP, Average
                for i, i_r in enumerate(i_radius):
                    new_arr[i_r] = new_arr[i_r] + lambdamap_v[i]
                    new_cnt[i_r] = new_cnt[i_r] + 1
            

        elif method == 'ellip':
            # npix = round(gal.nstar**(1/3) * 3) # to maintain roughly same pixel density. 
            # Or, constant physical scale
            import utils.sampling as smp
            import draw
            import matplotlib.pyplot as plt 
            from Cappellari import mge

            def find_nearest(array,value):
                idx = (np.abs(array-value)).argmin()
                return idx        
        
            npix_mge = npix# * 2 # twice more pixels to measure ellipses. No, let's keep it consistent
            region = smp.set_region(xc=0, yc=0, zc=0, radius = 0.5 * l_img)
            fig, ax = plt.subplots(1)
            
#            print("shape of data", data.shape)

            eps_arr = []
            mjr_arr = []
            pa_arr = []
            xpos_arr=[]
            ypos_arr=[]

            print("Reff (px)", reff / dx)
            for i in range(10):
                f = mge.find_galaxy.find_galaxy(self.mmap, quiet=True, plot=True,
                                                mask_shade=False,
                                                fraction=0.05 + 0.05 * i)
                mjr_arr.append(f.majoraxis)#f.majoraxis * 3.5 / npix * l_img
                eps_arr.append(f.eps)
                pa_arr.append(f.theta)
                xpos_arr.append(f.xmed)
                ypos_arr.append(f.ymed)
                #xpos_arr[i] = f.xpeak               #ypos_arr[i] = f.ypeak

            plt.savefig(galaxy_plot_dir + str(self.info.nout).zfill(3) \
                            + "_" + str(self.id) + "ellipticity.png")
            plt.close()


            mge_interpol = False
            
            if mge_interpol:
                def _interpol(x0, x1, reff, y0, y1 ):
                    return (y0 / (x1 - reff) + y1 / (reff - x0)) / (x1 - x0)

                dist_r = mjr_arr * np.sqrt(1 - eps_arr) * dx * 3.5 # Let's just assume..
#                print(reff, dist_r)
                ind = np.where(dist_r > reff)[0][0] -1 # Last element that is smaller than self.reff
                x0 = dist_r[ind]
                x1 = dist_r[ind+1]
                             
                eps = _interpol(x0, x1, reff, eps_arr[ind], eps_arr[ind+1])
                pa = _interpol(x0, x1, reff, pa_arr[ind], pa_arr[ind+1])
                sma = _interpol(x0, x1, reff, mjr_arr[ind], mjr_arr[ind+1])
                smi = sma * (1-eps)
                xcen = _interpol(x0, x1, reff, xpos_arr[ind], xpos_arr[ind+1])
                ycen = _interpol(x0, x1, reff, ypos_arr[ind], ypos_arr[ind+1])
                self.eps, self.pa, self.sma, self.smi = eps, pa, sma, smi

            else:
                i_reff = find_nearest(mjr_arr * np.sqrt(1 - np.square(eps_arr)), reff)
                #i_reff = 3
                print('i_reff', i_reff)
                self.eps = eps_arr[i_reff] # eps at 1 * R_half(=eff)
                self.pa  = pa_arr[i_reff]
                sma = reff / np.sqrt(1-self.eps) / dx
                smi = sma*(1-self.eps)
                self.mjr = sma
                self.mnr = self.mjr*(1-self.eps)
                xcen = xpos_arr[i_reff]
                ycen = ypos_arr[i_reff]
                #sma=mjr_arr[i_reff] * 0.5 * 3.5 # SEMI major axis, pixel unit

            pa_rad = -1 * self.pa / 180 * np.pi

            cos = np.cos(pa_rad)
            sin = np.sin(pa_rad)
            print("eps", self.eps)
            print("dx", dx)
            print("nx,ny", nx,ny)            
            print("sma, smi", sma, smi)
            print("xcen, ycen", xcen, ycen)
            
            dd = np.sqrt(((xNode-xcen)*cos + (yNode-ycen)*sin)**2/sma**2 + ((yNode-ycen)*cos - (xNode-xcen)*sin)**2/smi**2)

            i_close = np.where(dd < 1.0)[0]
            print("Reff = half light?1", sum(mmap[i_close])/ sum(mmap))


            dist1d = np.sqrt(np.square(xNode - xcen) + np.square(yNode - ycen))
            map2 = mmap.copy()
            ab = np.sqrt(sma*smi)
            for i in range(npoints):
                ind = np.where( (dd > i/reff) & (dd < (i+1)/reff))[0]
#                print(i*rr/npoints, len(ind))
                a = sum(mmap[ind] * dist1d[ind] * abs(vmap[ind]))
                if a != 0:
                    b = sum(mmap[ind] * dist1d[ind] 
                            * np.sqrt(np.square(vmap[ind]) + np.square(sigmap[ind])))
#                    ind2 = np.where(sigmap[ind] > 0)[0]
#                    b = sum(mmap[ind[ind2]] * dist1d[ind[ind2]] 
#                            * np.sqrt(vmap[ind[ind2]]**2 + sigmap[ind[ind2]]**2))
                    points[i] = a/b
            
            if voronoi:
                dd = np.sqrt( (xNode - 0.5*npix)**2 + (yNode - 0.5*npix)**2 )
                i_radius = np.fix(dd).astype(int)
                new_arr = np.zeros(np.fix(max(dd)) + 1)
                new_cnt = np.zeros(np.fix(max(dd)) + 1)
                # NGP, Average
                for i, i_r in enumerate(i_radius):
                    new_arr[i_r] = new_arr[i_r] + lambdamap_v[i]
                    new_cnt[i_r] = new_cnt[i_r] + 1
            
            # npix * npix map for plots
                for ibin in range(len(xNode)):
                    ind = np.where(binNum == ibin)[0]
                    vmap_org[ind] = vmap_v[ibin]
                    sigmap_org[ind] = sigmap_v[ibin]
            
                fig.suptitle("ID: {}    z: {:2f}".format(str(self.id).zfill(5),
                             self.info.zred))

# Stellar particle density map
                fig, axs = plt.subplots(2,2)
                axs[0,0].imshow(mmap_org.reshape(nx,ny),interpolation='nearest')
                axs[0,1].imshow(vmap_org.reshape(nx,ny),interpolation='nearest')
                axs[1,0].imshow(sigmap_org.reshape(nx,ny),interpolation='nearest')
                axs[1,1].plot(new_arr / new_cnt)
                axs[1,1].set_ylim(0,1)

            else:
                fig, axs = plt.subplots(2,2)
                cc = axs[0,0].imshow(dd.reshape(nx,ny))
                plt.colorbar(cc)
                axs[0,1].imshow(np.log10(mmap).reshape(nx,ny))
                axs[1,0].imshow(vmap.reshape(nx,ny))
                axs[1,1].imshow(sigmap.reshape(nx,ny))
                
            plt.savefig(galaxy_plot_dir + str(self.info.nout).zfill(3) \
                            + "_" + str(self.id) + "dd.png")
            plt.close()


        self.dist1d = dist1d
        self.lambda_arr = points
        # npix = ind_1Reff. 
        # ** npix 2 = npix * rscale
        # 0.5 * npix = Reff
        # 0.25 * npix = 0.5Reff
        self.lambda_r = np.average(self.lambda_arr[int(0.25 * npix) : int(0.25 * npix) + 1])
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
    
    def reorient(self, dest=[0., 0., 1], pop_nvec = ['star'], verbose=False):
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

        # New coordinates
        for target in ['star', 'dm', 'cell']:
            try:
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

            except:
                print("Error in reorient, couldn't update pos and vel")
                continue

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
        import numpy as np
                
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

        self.nvec = np.array([lx, ly, lz])/np.sqrt(lx**2 + ly**2 + lz**2)
#        self.cal_rotation_matrix(dest=dest)
        return self.nvec

    def _get_rotation_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        http://stackoverflow.com/a/6802723 
        https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
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
    
        vdiff = np.square(vx - self.vxc) + np.square(vy - self.vyc) + np.square(vz - self.vzc)
        npart = len(vdiff)
        
        if npart < 100:
            print("Too little particles within velocity window km/s")
            return False
    
        vdiff = np.zeros(npart, dtype=float)
        
        emin = np.argsort(vdiff)[:min(n, npart)]
        self.most_bound_particles = self.star['id'][emin]
        
    
    def cal_b2t(self, ptype='star', disk_criterion="Abadi2003",
                bound_only = True, hist=False):
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
            
        import numpy as np
        # Radius

        # Calculating boundness requires total mass inside a radius.
        # -> DM, Cell are also needed. 
        #part = getattr(self, ptype)
        
        m_g = gas_mass(self.cell, self.info)
        m_d = self.dm['m']# * info.msun
        
        m_all = np.concatenate((self.star['m'], m_g, m_d))

        r_s = np.sqrt(np.square(self.star['x']) 
                    + np.square(self.star['y']) 
                    + np.square(self.star['z']))
        r_g = np.sqrt(np.square(self.cell['x'])
                    + np.square(self.cell['y']) 
                    + np.square(self.cell['z']))
        r_d = np.sqrt(np.square(self.dm['x']) 
                    + np.square(self.dm['y']) 
                    + np.square(self.dm['z']))

        r_all = np.concatenate((r_s, r_g, r_d))
        i_sorted = np.argsort(r_all)
        m_enc = np.cumsum(m_all[i_sorted])

       
        # First nstar indices are stars.
        if bound_only:
            i_star = i_sorted[self.bound_ptcl] 
            x = self.star['x'][self.bound_ptcl]
            y = self.star['y'][self.bound_ptcl]
            vx = self.star['vx'][self.bound_ptcl]
            vy = self.star['vy'][self.bound_ptcl]
            m = self.star['m'][self.bound_ptcl]# * info.msun
        else:
            i_star = i_sorted[0:len(r_s)] 
            x = self.star['x']
            y = self.star['y']
            vx = self.star['vx']
            vy = self.star['vy']
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
            
        self.d2t = disk / self.mstar # or sum(self.star['m'][ind])
        if hist:
            return np.histogram(ellipticity, range=[-2,2], bins=20)



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
            
            
    def stellar_map(self, npix=400, proj='z'):
        """        
            returns stellar density map.
        """
        
        
         
    def plot_gal(self, npix=400, fn_save=None, ioff=True,
                 do9plots = False, **kwargs):
                
        import matplotlib.pyplot as plt
        from matplotlib import ticker
        from draw import pp
        from matplotlib.colors import LogNorm
#        import utils.prettyplot as ptt
        import numpy as np
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

        rgal_to_reff = self.Rgal_to_reff

        if do9plots:
            fig, axs = plt.subplots(3,3)
            fig.set_size_inches(12,9) # 9 plots
        else:
            fig, axs = plt.subplots(2,2)

        try:
            fig.suptitle("ID: {}    z: {:2f}".format(str(self.id).zfill(5), self.info.zred))
        except:
            fig.suptitle("ID: {}    z: not available".format(str(self.id).zfill(5)))

# Stellar particle density map
        ax = axs[0,0]
        """
            1st plot. - Stellar density map. 
            
        """        
        img = pp.part2den(self.star, self.info, npix=npix)

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

        plt.savefig(fn_save, dpi=200)
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
        outfile = hdf.File(base + str(self.id).zfill(6) + '_gal.hdf5',
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
        

    def load_gal(self, base='./'):
        """
        Load a galaxy from HDF5 file. 
        """
        pass#%%        
