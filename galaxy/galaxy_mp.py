# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 23:45:55 2015

Measure ICL fraction by merger channel


@author: hoseung

Inherites halo.HaloMeta. HaloMeta does not contain 
positional and kinematic information of the halo. This is good because
recalculated center of galaxy is likely to deviate from the center of the halo. 
"""

class Galaxy():

    def __init__(self, halo, radius_method='simple', info=None):
#        self.set_sim(sim)
        self.set_halo(halo)
        self.id = int(halo['id'])
        self.halo = halo
#    def set_sim(self, sim):
#        self.sim = sim
        self.info = info
        self.xc = 0
        self.yc = 0
        self.zc = 0
        self.reff = 0
        self.mstar = 0
        self.nstar = 0
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
    
    def mk_gal(self, star=None, dm=None, cell=None,
               save=False, rscale=0.3, verbose=False):
        import numpy as np

        print("Making a galaxy:", self.id)
        rscale_cen = 0.25
        xc_tmp = self.halo['x']
        yc_tmp = self.halo['y']
        zc_tmp = self.halo['z']
        rr_tmp = min([self.halo['rvir'] * rscale_cen, 0.0002]) # arbitrary! < 20kpc
        rr_tmp = max([rr_tmp, 0.000025])
        # When merger occurs, larger radius is likely to include 
        # companion galaxy resulting center to be in the middle of nowhere.
        # If you want larger galaxy, 
        # increase rgal_tmp instead. 
        # xx is easier to search for than x.

        if verbose:
            print("Searching particles insdie the radius..")
            print(rr_tmp * self.info.pboxsize * 1000, ['kpc'])
        
        xall = star['x']
        yall = star['y'] 
        zall = star['z'] 

        ind = np.where((xall - xc_tmp)**2 + (yall - yc_tmp)**2 
                            + (zall - zc_tmp)**2 < rr_tmp**2)[0]
        if verbose:
            print("there are {} particles inside the radius1".format(len(ind)))
        
        xc_tmp = np.median(xall[ind])
        yc_tmp = np.median(yall[ind])
        zc_tmp = np.median(zall[ind])
        if verbose:
            print("Now guessed center is :", xc_tmp, yc_tmp, zc_tmp)

        rr_tmp = min([self.halo['rvir'] * rscale, 0.00025]) # modified!
        if verbose:
            print("Search for particles around the guessed center, within...")
            print("radius2", rr_tmp * self.info.pboxsize * 1000, ['kpc'])

        ind = np.where((xall - xc_tmp)**2 + (yall - yc_tmp)**2 
                            + (zall - zc_tmp)**2 < rr_tmp**2)[0]
        if verbose:
            print("Number of particles:", len(ind))

        xx = xall[ind]
        yy = yall[ind] 
        zz = zall[ind] 
        mm = star['m'][ind] 

        if len(ind) < 100:
            
            self.star = False
            return False

# First, define the center w.r.t. stars.
        xc, yc, zc = self.get_center(xx, yy, zz, mm) # STARS!

#        self.physical_print([xc_tmp, yc_tmp, zc_tmp])

#        self.physical_print([xc, yc, zc])
        self.xc = xc
        self.yc = yc
        self.zc = zc

# Second, define the search radius
# Only 'simple' works, because 'star' is not defined yet.
        rgal_tmp = rr_tmp
        
        if dm is not None:
            idm = np.where( (dm["x"] - xc)**2 + 
                            (dm["y"] - yc)**2 + 
                            (dm["z"] - zc)**2 <= rgal_tmp**2)[0]
            ndm_tot = len(idm)

            self.dm = np.recarray(ndm_tot, dtype=dm.dtype)
            if 'id' in dm.dtype.names:
                self.dm['id'] = dm['id'][idm] 
            if 'm' in dm.dtype.names:
                self.dm['m'] = dm['m'][idm]
            if 'y' in dm.dtype.names:
                self.dm['x'] = dm['x'][idm] - self.xc
                self.dm['y'] = dm['y'][idm] - self.yc
                self.dm['z'] = dm['z'][idm] - self.zc
            if 'vy' in dm.dtype.names:
                self.dm['vx'] = dm['vx'][idm]# - self.vxc
                self.dm['vy'] = dm['vy'][idm]# - self.vyc
                self.dm['vz'] = dm['vz'][idm]# - self.vzc
            if verbose: print("DM data stored")
            
        if star is not None:
            istar = np.where((star["x"] - xc)**2 + 
                            (star["y"] - yc)**2 + 
                            (star["z"] - zc)**2 <= rgal_tmp**2)[0]
            nstar_tot = len(istar)
            if verbose: print("nstar tot:", nstar_tot)        
            if verbose: print("Store stellar particle")
    
            self.star = np.recarray(nstar_tot, dtype=star.dtype)
            if 'id' in star.dtype.names:
                self.star['id'] = star['id'][istar]
            if 'm' in star.dtype.names:
                self.star['m'] = star['m'][istar]
            if 'x' in star.dtype.names:
                self.star['x'] = star['x'][istar] - self.xc
            if 'y' in star.dtype.names:            
                self.star['y'] = star['y'][istar] - self.yc
            if 'z' in star.dtype.names:
                self.star['z'] = star['z'][istar] - self.zc
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

        # Or, hasattr(self, 'sink')
#        if 'sink' in part.pt:        
#            dtype = part._get_dtype("sink")
#            self.sink = np.recarray(nsink_tot * 2109, dtype=dtype)
#            isink = np.where( (part.star["x"] - xc)**2 + 
#                            (part.star["y"] - yc)**2 + 
#                            (part.star["z"] - zc)**2 <= rgal_tmp**2)[0]        
#            nsink_tot = len(isink)
        if cell is not None:
            icell = np.where((cell["x"] - xc)**2 + 
                (cell["y"] - yc)**2 + 
                (cell["z"] - zc)**2 <= rgal_tmp**2)[0]
            ncell_tot = len(icell)
            self.cell = np.recarray(ncell_tot, dtype=cell.dtype)
            self.cell['x'] = cell['x'][icell] - self.xc
            self.cell['y'] = cell['y'][icell] - self.yc
            self.cell['z'] = cell['z'][icell] - self.zc           
            self.cell['dx'] = cell['dx'][icell]
            self.cell['var0'] = cell['var0'][icell]
            self.cell['var2'] = cell['var2'][icell]
            self.cell['var3'] = cell['var3'][icell]
            self.cell['var4'] = cell['var4'][icell]
            self.cal_mgas()
            
        if verbose:
            print("Calculate Effective radius")
        self.reff = self.get_radius("eff")
        rgal_tmp = 4 * self.reff
        """
        print("Rgal = 4 * Reff = ", rgal_tmp * self.info.pboxsize * 1000)
        
            # Save sink particle as a BH, not cloud particles. 

        print('Simple check:')
        print("mean particle positions:", 
              np.mean(self.star['x']),
            np.mean(self.star['y']),
            np.mean(self.star['z']))


        print("median particle positions:", 
              np.median(self.star['x']),
            np.median(self.star['y']),
            np.median(self.star['z']))

        """
        self.nstar = nstar_tot
        self.mstar = sum(self.star['m'])
        # Now, get cov
        import utils.sampling as smp
        self.region = smp.set_region(xc=xc, yc=yc, zc=zc, radius = rgal_tmp)
        
#        print("before get_cov", time.time() - t0)
        self.get_cov(center_only=True)
#        print("after get_cov", time.time() - t0)
        return True

    def get_radius(self, method):
        import numpy as np
        if method == "simple":
            return self.halo['rvir']
        elif method == "eff":
            # requires galaxy center to be determined. 
            dd = self.star['x']**2 + self.star['y']**2 + self.star['z']**2 
#            npart = len(dd)
            dmax = max(dd)
            r_again = dmax
            m_annuli=[]
            nbin=8
            for i in range(nbin):
                m_annuli.append(sum(
                self.star['m'][(dd > i * dmax/nbin) * (dd < (i+1)*dmax/nbin)]))

            for i in range(nbin):
                if m_annuli[i] > m_annuli[i-1]:
                    r_again = dmax * (i+1)/nbin
                    break
            i_again = dd < r_again
            dsort = np.argsort(dd[i_again])
            m = self.star['m'][i_again]
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
        IDL version was to iteratively search for mass weighted mean position of
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

    def get_cov2(self):
        self.vxc = self.halo['vx']
        self.vyc = self.halo['vy']
        self.vzc = self.halo['vz']
        
        self.star['vx'] -= self.vxc
        self.star['vy'] -= self.vyc
        self.star['vz'] -= self.vzc


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
        
        # Add adaptivity
        if method == "mostbound":
            while npart > 80000:
                r_frac *= 0.7
                vind = np.where((vx - vxt)**2 + 
                    (vy - vyt)**2 +
                    (vz - vzt)**2 < r_frac * vrange**2)[0]
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
            multithreaded=False
            if multithreaded:
                from multiprocessing import Pool
                from itertools import repeat
                """ Why no performance improvement at all..???"""
        
                pool=Pool(processes=4)
                e = pool.starmap(self.e_vdiff,
                                 zip(range(1, len(vxt)-2),
                                    repeat(vxt),
                                    repeat(vyt),
                                    repeat(vzt)))
            elif not multithreaded:
                for iv in range(1, len(vxt)-2):
                    vdiff=np.sqrt((vxt[0:-1:5] - vxt[iv])**2 + 
                                  (vyt[0:-1:5] - vyt[iv])**2 +
                                  (vzt[0:-1:5] - vzt[iv])**2)
                    e[iv] = sum(-1./vdiff[vdiff != 0])
    
# N parts that are most bound to <50000 'central' particle group.
            emin = np.argsort(e)[:5000]
            
            m = self.star['m'][vind]
            self.vxc = np.average(vxt[emin], weights=m[emin]) # top 100 
            self.vyc = np.average(vyt[emin], weights=m[emin])
            self.vzc = np.average(vzt[emin], weights=m[emin])
        elif method == 'com':
            self.vxc, self.vyc, self.vzc = self.get_center(
                                vx, vy, vz, self.star['m'], tol=5, niter=8)
        
        self.star['vx'] -= self.vxc
        self.star['vy'] -= self.vyc
        self.star['vz'] -= self.vzc
    
    def cal_mgas(self):
        self.mgas = 0.0
        pass
    
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
    
        
    def cal_lambda_r(self, npix=30, rscale=3.0, method=2):
        import numpy as np

        # Only within effective radius. (or whatever radius.)
        # already centered.
        self.rscale_lambda = rscale
        rr = self.reff * rscale
        npix2 = round(npix * rscale)

        ind_ok = np.where((abs(self.star["x"]) <= rr) &
                          (abs(self.star["y"]) <= rr) &
                          (abs(self.star["z"]) <= rr) )[0]
# Some margin make mmap and vmap look better.
# But only inside Lambda_R for R < Reff is returned. 

        print("{} particles inside {}kpc".format(len(ind_ok), 
                    rr * self.info.pboxsize * 1000))

        x = self.star['x'][ind_ok]
        y = self.star['y'][ind_ok]
        m = self.star['m'][ind_ok]
        vz = self.star['vz'][ind_ok]
        
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
        
        dx = (2 * rr * self.info.pboxsize * 1000) / npix2
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

        points = np.zeros(0.5 * npix) # 1.5 Reff
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
        self.lambda_r = 0.5*(self.lambda_arr[0.5*len(self.lambda_arr) -1] 
                                + self.lambda_arr[len(self.lambda_arr) -1])
        if verbose: print("lambda_arr done")
    
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
        
    def bound_particles(self, ptype='star'):
        """
        Returns indices of bound particles. 
        By default for stellar particles.
        """       
    
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
            import matplotlib.pyplot as plt
            plt.hist(ellipticity, bins=20, range=[-2,2])
            return np.histogram(ellipticity, range=[-2,2], bins=20)


    def cal_trivia(self):
        self.mstar = sum(self.star['m'])
        self.vcen = self.get_vcen()
#        self.metal = 
#        self.reff = 12333

    def save_gal_pickle(self):
        import pickle
        with open(str(self.id).zfill(6) + 'gal.pickle', 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
         
    def plot_gal(self, npix=400, wdir='./', fname=None, **kwargs):
        import matplotlib.pyplot as plt
        from matplotlib import ticker
        from draw import pp
        from matplotlib.colors import LogNorm
#        import utils.prettyplot as ptt
        import numpy as np
        
        def fmt(x, pos):
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)

        do9plots = False        
        if do9plots:
            fig, axs = plt.subplots(3,3)
            fig.set_size_inches(12,9) # 9 plots
        else:
            fig, axs = plt.subplots(2,2)
        
        fig.suptitle("ID: {}    z: {:2f}".format(str(self.id).zfill(5), self.info.zred))

# Stellar particle density map
        ax = axs[0,0]
        rgal = self.reff * self.rscale_lambda

        img = pp.part2den(self.star, self.info, npix=npix)
#        region=self.region, offset=[self.xc, self.yc, self.zc])
#        print("img max:", max(img.data))
        im = ax.imshow(img.data, cmap=plt.get_cmap('brg'),
                       norm=LogNorm(vmin=1e6))
        cb = plt.colorbar(im, ax=ax)

        ax.set_xlabel("position ["+ r'$R/R_{eff}$'+"]")
        ax.set_xticks([0, 100, 200, 300, 400])
        ax.set_xticklabels(["-4", "-2", "0", "2", "4"])
        ax.set_ylabel("position [kpc]")
        ax.set_yticks([0, 100, 200, 300, 400])
        yticks = ["{:.2f}".format(y * self.info.pboxsize * 1000) \
                    for y in np.linspace(-rgal, rgal, num=5)]
        ax.set_yticklabels(yticks)

# Lambda_r plot
        if self.lambda_arr is not None:
            ax = axs[0,1]
            ax.plot(self.lambda_arr) # ~ 1 * Reff
            ax.set_title(r"$\lambda _{R}$") 
            ax.text(0.5 * len(self.lambda_arr), 0.8, "{:.2e}".format(self.mstar) + r"$M_{\odot}$")
            ax.set_ylim(bottom=0, top=1)
            ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
            ax.set_xticks(np.linspace(0, len(self.lambda_arr), num=3))
            ax.set_xticklabels(["0", "0.5", "1"])


        """
        im = ax.imshow(self.mmap, cmap=plt.get_cmap('brg'),
                       norm=LogNorm(vmin=1e6))

        ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
        ax.set_xticks(np.linspace(0, len(self.mmap), num=3))
        ax.set_xticklabels(["-1.5", "0", "1.5"])
        ax.set_ylabel("[kpc]")
        ax.set_yticks(np.linspace(0, len(self.mmap), num=5))
        yticks = ["{:.2f}".format(y * self.info.pboxsize * 1000) \
                    for y in np.linspace(-self.reff, self.reff, num=5)]
        ax.set_yticklabels(yticks)        

        cb = plt.colorbar(im,
                          ax=ax,
#                          format=ticker.FuncFormatter(fmt),
                          label=r"$M_{\odot} kpc^{-2}$")
#        tick_locator = ticker.MaxNLocator(nbins=3)
#        cb.locator = tick_locator
#        cb.update_ticks()
        ax.set_title("Stellar density map")
        """        
# sigma map
        if hasattr(self, "sigmap"):
            ax = axs[1,0]        
            im = ax.imshow(self.sigmap, vmin=0, vmax = 200, cmap=plt.get_cmap('brg'))
            cb = plt.colorbar(im, ax=ax, label=r'$km s^{-1}$') # needs a mappable object.
            tick_locator = ticker.MaxNLocator(nbins=3)
            cb.locator = tick_locator
            cb.update_ticks()
            ax.set_title(r"$\sigma$ map")
            ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
            ax.set_xticks(np.linspace(0, len(self.mmap), num=3))
            ax.set_xticklabels(["-1", "0", "1"])
            ax.set_ylabel("[kpc]")
            ax.set_yticks(np.linspace(0, len(self.mmap), num=5))
            yticks = ["{:.2f}".format(y * self.info.pboxsize * 1000) \
                        for y in np.linspace(-self.reff, self.reff, num=5)]
            ax.set_yticklabels(yticks)        
#        ax.locator_params(tight=True, nbins=5)

# velocity map
        if self.vmap is not None:
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
            yticks = ["{:.2f}".format(y * self.info.pboxsize * 1000) \
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
        if fname is None:
            fname = "galaxy_plot" + str(self.id).zfill(5) + ".png"
        
        plt.savefig(wdir + fname, dpi=100)
        plt.close()
        
        
    def save_gal(self, base='./'):
        import h5py as hdf
        # Save data into a hdf5 file
        filename = base + str(self.id).zfill(6) + 'gal.hdf5'
        
        # Create and open the file with the given name
        with hdf.File(filename, 'w', libver='latest') as outfile:
        
        # Make a data set called galaxy
#        outfile.create_dataset("galaxy")
        
        # Create a group for stroing data
            grp_gal = outfile.create_group("galaxy")
    
            # Store data under /galaxy with direct assignment        
            grp_gal.create_dataset("star/x", data=self.star['x'])
            grp_gal.create_dataset("star/y", data=self.star['y'])
            grp_gal.create_dataset("star/z", data=self.star['z'])
            grp_gal.create_dataset("star/m", data=self.star['m'])
            grp_gal.create_dataset("star/id", data=self.star['id'])
            grp_gal.create_dataset("star/vx", data=self.star['vx'])
            grp_gal.create_dataset("star/vy", data=self.star['vy'])
            grp_gal.create_dataset("star/vz", data=self.star['vz'])
            if "time" in self.pq:
                grp_gal.create_dataset("star/time", data=self.star['time'])
            if "metal" in self.pq:
                grp_gal.create_dataset("star/metal", data=self.star['metal'])        
            
            if "dm" in self.pq:
                grp_gal.create_dataset("dm/x", data=self.dm['x'])
                grp_gal.create_dataset("dm/y", data=self.dm['y'])
                grp_gal.create_dataset("dm/z", data=self.dm['z'])
                grp_gal.create_dataset("dm/m", data=self.dm['m'])
                grp_gal.create_dataset("dm/id", data=self.dm['id'])
                grp_gal.create_dataset("dm/vx", data=self.dm['vx'])
                grp_gal.create_dataset("dm/vy", data=self.dm['vy'])
                grp_gal.create_dataset("dm/vz", data=self.dm['vz'])
        """
        grp_gal.create_dataset("cell/x", data=self.cell['x'])
        grp_gal.create_dataset("cell/y", data=self.cell['y'])
        grp_gal.create_dataset("cell/z", data=self.cell['z'])
        grp_gal.create_dataset("cell/vx", data=self.cell['var1'])
        grp_gal.create_dataset("cell/vy", data=self.cell['var2'])
        grp_gal.create_dataset("cell/vz", data=self.cell['var3'])
        grp_gal.create_dataset("cell/rho", data=self.cell['var0'])
        grp_gal.create_dataset("cell/temp", data=self.cell['var4'])
        grp_gal.create_dataset("cell/metal", data=self.cell['var5'])
        # Check the units, especially for temperature and metallicity
        """
        
    def load_gal(self, base='./'):
        """
        Load a galaxy from HDF5 file. 
        """
        pass
#%%        
