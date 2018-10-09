import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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


def gen_vmap_sigmap(self,
                    npix_per_reff=5,
                    rscale=3.0,
                    n_pseudo=1,
                    voronoi=None,
                    verbose=False,
                    plot_map=False,
                    weight="m"):
    """
    generates mass, velocity, sigma map from stellar particles.
    npix = (npix_per_reff * 2 * (rsacle + 1))
    and mmap, vmap, sigmap are npix x npix 2d arrays.
    So, if rscale = 3 is given, calculate inside 3Reff but plot a map of 4Reff.


    Parameters
    ----------
    npix_per_reff: int
        number of bins/pixels per 1reff. default = 5
    rscale : float
        size of galaxy = rscale * reff
    n_pseudo : int
        split each particle into n_pseudo pseudo particles. default = 1
    voronoi :
        parameter set of voronoi tesselation as a dict. None by default.
    weight : field name
        Name of the field to give wieght. for example:  "m", "flux_g", "flux_r",...
    verbose :
        default = False


    Todo
    1. if n_pseudo > 1, sig = 0.3kpc should scale with aexp.
    2. Separate lambda calculation with voronoi tesselation.

    Note
    ----
    As an analogy to np.histogram2d, there is a convenience function binned_statistic_2d.
    For example, you can calculate std of each pixel as
    >>> sig_map = binned_statistic_2d(gg.star["x"], gg.star["y"], gg.star["vz"], statistic=np.std, bins=40)

    However, this convenience function is slower than this gen_vmap_sigmap() function (all packages updated at 2017.07.01).
    ( > 18s Vs 10s for 60*1e5 particles, > 0.6s Vs 0.5s for 1e5 particles)
    So, don't bother looking for other implementations ;)

    """
    # already centered.
    reff = self.meta.reff
    # reff restriction must be given at earlier stage, radial_profile_cut()
    r_img_kpc = reff * (rscale+1) # in kpc
    # Some margin makes mmap and vmap look better.
    # If rscale = 3 is given, calculate inside 3Reff,
    # but plot a map of 4Reff.

    dx = reff / npix_per_reff # kpc/pixel
    npix = round(npix_per_reff * 2 * (rscale+1))
    nx, ny = npix, npix
    # to suppress contamination from tidal tail,
    # give a cylindrical cut, not spherical cut.
    # Cappellari 2002 assumes Cylindrical velocity ellipsoid.
    # particles inside 4Reff.
    ind = (np.square(self.star['x']) + \
           np.square(self.star['y'])) < np.square(reff * (rscale + 1))

    n_frac = sum(ind)/self.meta.nstar*100.0
    if verbose :
        print("{:.2f}% of stellar particles selected".format(n_frac))
    if n_frac < 10:
        print("Too few stars are selected...")
        print("min max x", min(self.star['x']), max(self.star['x']))
        print("min max y", min(self.star['y']), max(self.star['y']))
        print("min max z", min(self.star['z']), max(self.star['z']))
        print("# star", len(ind))
        return [-1,-1,-1], [-1,-1,-1]
    # 100% means that the galaxy radius is smaller than 4Reff.
    # Considering ellipticity, even 4Reff may not be enough to derive.

    if n_pseudo > 1:
        # sig in kpc unit. up to 1M particles
        # sig = 0.3kpc from Naab 2014.
        # Todo
        # sig = 0.3kpc should scale with aexp.
        n_pseudo = max([round(1e6/self.meta.nstar), n_pseudo])
        xstars, ystars, mm, vz = self._pseudo_particles(self.star['x'][ind],
                                                  self.star['y'][ind],
                                                  self.star[weight][ind],
                                                  self.star['vz'][ind],
                                                  sig=0.3,
                                                  n_times=n_pseudo)
    else:
        xstars = self.star['x'][ind]
        ystars = self.star['y'][ind]
        vz = self.star['vz'][ind]
        if self.star[weight][0] ==0:
            print("[galaxymodule.rotation_parameter] warning... \n"+
                  "the weight {} is zero! Fall back to mass".format(weight))
            mm = self.star["m"][ind]
        else:
            mm = self.star[weight][ind]

    if verbose: print(("\n" "Calculating rotation parameter using {} particles "
    "inside {:.3f}kpc, or {}Reff".format(len(ind), r_img_kpc, rscale + 1)))

    # using NGP charge assignment
    # fix center explicitly.
    # 0.5 * (min + max) != center
    xstars = (xstars + r_img_kpc) / r_img_kpc * 0.5 * nx # 0 < xstarts < nx
    ystars = (ystars + r_img_kpc) / r_img_kpc * 0.5 * ny
    # because of Gaussian smoothing, pseudo particles can go out of cic region.
    # But don't worry, np.clip is ready.

    # NGP assignment
    ngx = np.clip(np.fix(xstars), 0, nx-1)
    ngy = np.clip(np.fix(ystars), 0, ny-1)
    indices = (ngx + ngy * nx).astype(np.int32)

    # Mass map
    mmap = np.zeros(nx * ny, dtype=float) # should cover 4Reff.
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

    if voronoi is not None:
        noise_map = np.sqrt(count_map)
        noise_map[noise_map < 1] = 1 # minimum noise for empty pixeles

        xpos_regular = np.repeat(np.arange(nx),ny)
        ypos_regular = np.tile(np.arange(ny),nx)
        from Cappellari.voronoi.voronoi_2d_binning import voronoi_2d_binning
        """
        This function accepts only data on uniform grid...?
        => Yes.
        Provide pixelsize to speed up the calculation.
        """
        print("Inside voronoi, plot =", voronoi["plot"])
        binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = \
            voronoi_2d_binning(xpos_regular, ypos_regular, count_map,
                             noise_map, pixelsize=1, targetSN=voronoi["targetSN"],
                             plot=voronoi["plot"], quiet=voronoi["quiet"])

        self.xNode = xNode
        self.yNode = yNode
        self.binNum = binNum
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

            # mass-weighted sigma
            sigmap_v[ibin] = self.weighted_std(vz[i_part], weights=mm[i_part])

            # update original map too.
            #vmap[ind] = vmap_v[ibin]
            #sigmap[ind] = sigmap_v[ibin]
        self.mmap_v = mmap_v
        self.vmap_v = vmap_v
        self.sigmap_v = sigmap_v

        # Always make 2D maps.
        vmap_plot = np.empty_like(count_map).ravel()
        sigmap_plot = np.empty_like(count_map).ravel()
        # npix * npix map for plots
        for ibin in range(len(xNode)):
            ind = np.where(binNum == ibin)[0]
            vmap_plot[ind] = vmap_v[ibin]
            sigmap_plot[ind] = sigmap_v[ibin]
        # Stellar particle density map
        self.mmap = mmap.reshape(nx,ny)
        self.vmap = vmap_plot.reshape(nx,ny)
        self.sigmap = sigmap_plot.reshape(nx,ny)

    else:
        # Velocity and dispersion map
        vmap = np.zeros(nx * ny, dtype=float)
        sigmap=np.zeros(nx * ny, dtype=float)

        # No Voronoi tessellation.
        for i in range(nx * ny):
            ind = np.where(indices == i)[0]
            if len(ind) > 0:
                # mass-weighted sigma
                sigmap[i] = self.weighted_std(vz[ind], mm[ind])
                # mass-weighted velocity
                vmap[i] = np.average(vz[ind], weights=mm[ind])
            else:
                sigmap[i] = 0
                vmap[i] = 0
        #self.xNode = np.tile(np.arange(nx),ny) # x = tile? or repeat?
        #self.yNode = np.repeat(np.arange(ny),nx)
        self.sigmap = sigmap.reshape(nx, ny)
        self.vmap = vmap.reshape(nx,ny)

    if plot_map:
        fig, axs = plt.subplots(2,2)
        axs[0,0].imshow(self.mmap,interpolation='nearest', norm=LogNorm())
        axs[0,1].imshow(self.vmap,interpolation='nearest')
        axs[1,0].imshow(self.sigmap,interpolation='nearest')
        axs[1,1].set_ylim(0,1)
        plt.savefig("vsigmaps_{}.png".format(self.meta.id))
        plt.close()


def get_mge_out(f, frac, npix_per_reff, name):
    sma = npix_per_reff
    pa_rad = -1*f.theta/180*np.pi
    return dict({"name":name,
                 "frac":frac,
                 "eps":f.eps,
                 "sma":sma,
                 "smi":sma * (1-f.eps),
                 "pa":f.theta,
                 "pa_rad":pa_rad,
                 "cos":np.cos(pa_rad),
                 "sin":np.sin(pa_rad),
                 "xcen":f.xmed,
                 "ycen":f.ymed})


def _measure_lambda(mge_par,
                    mmap, vmap, sigmap,
                    xNode, yNode,
                    npix_per_reff,
                    rscale,
                    voronoi=False,
                    do_plot=False,
                    verbose=False,
                    ):

    xcen=mge_par["xcen"]
    ycen=mge_par["ycen"]
    cos=mge_par["cos"]
    sin=mge_par["sin"]
    sma=mge_par["sma"]
    smi=mge_par["smi"]

    dd = np.sqrt(((xNode-xcen)*cos + (yNode-ycen)*sin)**2/sma**2 + \
                 ((yNode-ycen)*cos - (xNode-xcen)*sin)**2/smi**2) * \
                 npix_per_reff

    # lambda calculaed over '3'Reff.
    npix = round(npix_per_reff * rscale)
    points = np.zeros(npix)

    if verbose: print("Reff = half light?1", sum(mmap[dd < 1.0])/ sum(mmap))
    dist1d = np.sqrt(np.square(xNode - xcen) + np.square(yNode - ycen))

    for i in range(len(points)):
        #ind = np.where( (dd > i) & (dd < (i+1)))[0]
        ind = np.where(dd < (i+1))[0]

        if len(ind) > 0:
            a = np.sum(mmap[ind] * dist1d[ind] * abs(vmap[ind]))
            if a > 0:
                ind2 = np.where(sigmap[ind] > 0)[0]
                b = np.sum(mmap[ind[ind2]] * dist1d[ind[ind2]]
                        * np.sqrt(vmap[ind[ind2]]**2 + sigmap[ind[ind2]]**2))

                points[i] = a/b

    return points


def interpol(x, x0, arrays):
    ind = max([1,np.argmax(x > x0)])
    xl, xr = x[ind -1], x[ind]
    fl = (x0 - xl)/(xr-xl)
    fr = (xr - x0)/(xr-xl)

    return [y[ind -1]*fr + y[ind]*fl for y in arrays]


def cal_lambda_r_eps(self,
                     method='ellip',
                     verbose=False,
                     galaxy_plot_dir='./',
                     save_result = True,
                     recenter_v=True,
                     iterate_mge = False):
    """
    Calculate galaxy rotation parameter $\lambda$ following the definition of
    Emsellem 2007.

    This must be called after preparing velocity map and sigma map of the galaxy,
    which are generated by gen_vmap_sigmap(). It also relies on parameters used
    in gen_vmap_sigmap() method.


    """

    npix_per_reff=self.params.vmap_sigmap["npix_per_reff"]
    rscale=self.params.vmap_sigmap["rscale"]
    voronoi=self.params.vmap_sigmap["voronoi"]

    import matplotlib.pyplot as plt
    is_voronoi = lambda x : x is not None

    # 2. calculate profile over radial bins.
    # iterate_mge : run MGE with various light fraction
    #               to find out half light radius.
    # mge_interpol : from iterative measruements of MGE,
    #                interpolate quantities at the half light radius.
    mmap = np.ravel(self.mmap)

    nx = round(npix_per_reff * 2 * (rscale + 1))
    ny = nx

    xpos = np.tile(np.arange(nx),ny)
    ypos = np.repeat(np.arange(ny),nx)

    #xNode = self.xNode
    #yNode = self.yNode
    reff = self.meta.reff
    self.npix_per_reff = npix_per_reff

    mmap_tot = sum(mmap)
    dist = np.sqrt(np.square(xpos- nx/2) \
                 + np.square(ypos- ny/2))

    if method == 'ellip':
        from Cappellari import mge
        if iterate_mge:
            eps_arr = []
            mjr_arr = []
            pa_arr  = []
            xpos_arr= []
            ypos_arr= []
            f_light_arr=[]
            for i in range(6):
                # mmap = 1D, self.mmap = mmap.reshap(nx,ny)
                f = mge.find_galaxy.find_galaxy(self.mmap, quiet=True, plot=False,
                                                mask_shade=True,
                                                fraction=0.04*(i+1)**1.5)
                mjr_arr.append(f.majoraxis)#f.majoraxis * 3.5 / npix * l_img
                eps_arr.append(f.eps)
                pa_arr.append(f.theta)
                xpos_arr.append(f.xmed)
                ypos_arr.append(f.ymed)
                sma = f.majoraxis# * 3.5
                smi = sma * (1-f.eps)

                pa_rad = -1 * f.theta / 180 * np.pi
                cos = np.cos(pa_rad)
                sin = np.sin(pa_rad)

                f_light_arr.append(sum(mmap[((xNode-f.xmed)*cos + (yNode-f.ymed)*sin)**2/sma**2 + \
                                        ((yNode-f.ymed)*cos - (xNode-f.xmed)*sin)**2/smi**2 \
                                        < 3.5])/ mmap_tot)
            print("Reff = half light?, mjr, epss", f_light_arr[i], f.majoraxis, f.eps)

# Determine eps, pa, sma, xcen, ycen
            if mge_interpol:

                self.meta.eps, self.meta.pa, self.meta.sma, self.meta.xcen, self.meta.ycen = \
                    interpol(np.array(f_light_arr), 0.5, (eps_arr, pa_arr, mjr_arr, xpos_arr, ypos_arr))
                print('eps array', eps_arr)
                print("interpolated eps, pa, sma, xcen, ycen", \
                      self.meta.eps, self.meta.pa, self.meta.sma, self.meta.xcen, self.meta.ycen)
                self.meta.smi = self.meta.sma * (1 - self.meta.eps)

            else:
                i_reff = np.argmax(np.array(f_light_arr) > 0.5) -1# first element > 0.5 -1
                self.meta.eps = eps_arr[i_reff] # eps at 1 * R_half(=eff)
                self.meta.pa  = pa_arr[i_reff]
                sma = reff# / np.sqrt(1-self.meta.eps) / dx
                # sma becomes too large with large e, (when e = 0.9, sma = 10 * reff).
                smi = sma*(1-self.meta.eps)
                self.meta.sma = mjr_arr[i_reff] * 3.5
                self.meta.smi = self.meta.sma*(1-self.meta.eps)
                #xcen = xpos_arr[i_reff]
                #ycen = ypos_arr[i_reff]
                #sma=mjr_arr[i_reff] * 0.5 * 3.5 # SEMI major axis, pixel unit
        else:
            # MGE in one go.
            # mmap = 1D, self.mmap = mmap.reshap(nx,ny)
            # No iteration, 50% light in one shot.
            dsort = np.argsort(dist)
            # level = pixel value at cumsum = 50%.
            i_reff = np.argmax(np.cumsum(mmap[dsort]) > 0.5*mmap_tot)
            frac1 = i_reff/ len(mmap)

            d_05reff = np.argmax(dsort > 0.5*dsort[i_reff]) # dsort[d_05reff] = 0.5Reff
            frac05 = d_05reff/len(mmap)

            frac15 = np.pi / self.meta.rscale_lambda**2 * (10/(2*reff))**2
            # fraction of 15kpc in the mmap.
            frac15 = min([0.99, frac15])

            fracs = [frac1, frac05, frac15]
            names = ["1Reff", "0.5Reff", "15kpc"]
            self.meta.mge_result_list=[]
            self.meta.lambda_result_list=[]
            self.meta.lambda_r=[]
            for frac, name in zip(fracs, names):
                if verbose: print("frac", frac)
                f = mge.find_galaxy.find_galaxy(self.mmap, quiet=True, plot=False,
                                            mask_shade=False,
                                            fraction=frac)
                mge_now = get_mge_out(f, frac, npix_per_reff, name)
                # recenter velocity
                if recenter_v:
                    xcen=mge_now["xcen"]
                    ycen=mge_now["ycen"]
                    if is_voronoi(voronoi):
                        dist1d = np.sqrt(np.square(self.xNode - xcen) + np.square(self.yNode - ycen))
                        i_very_cen = np.where(dist1d < 1.5)[0] # within 0.2 Reff
                        #print("{:1f} % are very close to the center".format(100*len(i_very_cen)/len(dist1d)))
                        new_vc = np.mean(self.vmap_v[i_very_cen])
                        self.vmap_v = self.vmap_v - new_vc
                    else:
                        dist = np.sqrt(np.square(xpos- xcen) \
                                       + np.square(ypos- ycen))
                        i_very_cen = np.where(dist < 1.5)[0] # within 0.2 Reff
                        new_vc = np.mean(self.vmap.ravel()[i_very_cen])
                        self.vmap = self.vmap - new_vc

                self.meta.mge_result_list.append(mge_now)

                if is_voronoi(voronoi):
                    larr = _measure_lambda(mge_now,
                                       self.mmap_v, self.vmap_v, self.sigmap_v,
                                       self.xNode, self.yNode,
                                       npix_per_reff,
                                       rscale,
                                       voronoi=True,
                                       verbose=verbose)
                else:
                    larr = _measure_lambda(mge_now,
                                       mmap.ravel(), self.vmap.ravel(), self.sigmap.ravel(),
                                       xpos, ypos,
                                       npix_per_reff,
                                       rscale,
                                       voronoi=False,
                                       verbose=verbose)
                self.meta.lambda_result_list.append(larr)
                self.meta.lambda_r.append(np.average(larr[npix_per_reff - 1 : npix_per_reff + 2]))
