# -*- coding: utf-8 -*-
"""
Test bed for galaxy related routines in serial!

Created on Sun Aug 16 16:49:47 2015

@author: hoseung
"""

import numpy as np
import matplotlib.pyplot as plt
import load
import time

def prepare_data(base='./', nout=187, npix=800, lmax=19, idgal=1197):
    import tree.halomodule as hmo
    import utils.sampling as smp    
#    r_cluster_scale = 2.0
    r_gal_scale = 1.5

    ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]
    
    info = load.info.Info()
    info.setup(nout=nout, base=wdir)
    info.load()
    hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info)
    hh.load()
#    halo = hh.data
    
    # List of galaxies
#    i_center = np.where(halo['np'] == max(halo['np']))[0]
#    i_satellites = smp.extract_halos_within(hh.data,
#                    i_center, scale=r_cluster_scale, Mcut=1e5)
#    print("Total {0} halos \n{1} halos are selected".format(
#          len(i_satellites),sum(i_satellites)))
#    print(i_satellites[195:197])

    # Choose one halo
    # 5427(1198)
    h = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info)
    h.derive_from(hh, idgal)
    print(h.data['x'])
    region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z
                                  , radius = h.data.rvir * r_gal_scale)
    print(region)
    s = load.sim.Sim()
    s.setup(nout, wdir)    
    s.set_ranges(region["ranges"])
    s.show_cpus()
    s.add_part(ptypes)    
    s.add_hydro()
    t0 = time.time()
    s.part.load(fortran=True)
    t1 = time.time()
    s.hydro.amr2cell(lmax=19)
    t2 = time.time()
    print("Loading particle took {}, \n and loading hydro took {}".format(t1-t0, t2-t1))    

    return h, info, s.part.star, s.part.dm, s.hydro.cell

#%%------------------------------------------------------------------------
def _display_pixels(x, y, counts, pixelSize=None):
    """
    Display pixels at coordinates (x, y) coloured with "counts".
    This routine is fast but not fully general as it assumes the spaxels
    are on a regular grid. This needs not be the case for Voronoi binning.

    """
    from scipy.spatial import distance
    import numpy as np
    if pixelSize is None:
        pixelSize = np.min(distance.pdist(np.column_stack([x, y])))
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



def plot_halo_region(halodata, info, star, cell, npix=400, save_dir='./'):
    import draw

    import utils.sampling as smp
    
    extent = (0, npix, 0, npix)
    region = smp.set_region(xc = halodata['x'],
                            yc = halodata['y'],
                            zc = halodata['z'],
                            radius = halodata['rvir'])

    star_map = draw.pp.den2d(star['x'],star['y'],star['z'],star['m'], npix,
                             region=region, cic=True, norm_integer=False)
    if star_map is not False:
        ls = np.zeros((npix,npix))
        ii = star_map > 0
        ls[ii] = np.log10(star_map[ii]) # Stellar map HAS empty pixels.
        ls[star_map <= 0] = np.floor(ls.min())
        im1 = plt.imshow(ls, cmap="gray", interpolation='nearest', extent=extent)       
    
    # One of two should be transposed.
    # But which one?
    gas_map = draw.pp.pp_cell(cell, npix, info, region=region, verbose=False)
    im2 = plt.imshow(np.transpose(np.log10(gas_map)), cmap="CMRmap", alpha=.5, interpolation='bilinear', extent=extent)

    rgal = region['radius'] * info.pboxsize * 1000

    ax = plt.gca()
    ax.set_xlabel("position [kpc]")
    ax.set_xticks(np.linspace(0,npix,5))
    xticks = ["{:.2f}".format(x) \
                for x in np.linspace(-rgal, rgal, num=5)]
    ax.set_xticklabels(xticks)
    ax.set_ylabel("position [kpc]")
    ax.set_yticks(np.linspace(0,npix,5))
    yticks = ["{:.2f}".format(y) \
                for y in np.linspace(-rgal, rgal, num=5)]
    ax.set_yticklabels(yticks)
    
    plt.savefig(save_dir+"2dmap_"+str(halodata['id']).zfill(5)+'.png', dpi=144)
    plt.close()

def mk_gal(halodata, info,
           save=False, rscale=0.3, verbose=False, save_dir='./',
           rscale_lambda=2.0, npix=400):
    """
    Direct plot,
    Create galaxy, 
    Calculate lambda_r (using Cappellari 2003)
    Draw ed map of galaxy.

    """
    from galaxymodule import galaxy

    #Create galaxy ---------------------------------------------------------
    gal = galaxy.Galaxy(halodata, radius_method='simple', info=info)
    is_gal = gal.mk_gal(star=star, dm=dm, cell=cell,
               rscale=rscale, verbose=verbose)
    #-----------------------------------------------------------------------    
    print(gal.id, "IS_GAL",is_gal)
    if not is_gal:
        print(gal.id, " Not a good galaxy")
        # Save to catalog ----------------------------------------------------

    return gal
    

def plot_lambda(catalog, i_early, i_late, i_bad, out_dir='./'):
    import matplotlib.pyplot as plt
    plt.ioff()
    f = plt.figure()
    ax = f.add_subplot(111)
    #for i, val in enumerate(lambdar_arr):
    for i in i_early:
        a = np.asarray(catalog['lambda_arr'][i])
        ax.plot(a, 'r-', alpha=0.5) # Red = Early
    for i in i_late:
        ax.plot(catalog['lambda_arr'][i], 'b-', alpha=0.3) # Red = Early
    
    #plt.xlabel() # in the unit of Reff
    ax.set_title(r"$\lambda _{R}$") 
    ax.set_ylabel(r"$\lambda _{R}$") 
    ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
    ax.set_xlim(right=9)
    ax.set_xticks([0, 4.5, 9])
    ax.set_xticklabels(["0", "0.5", "1"])
    plt.savefig(out_dir + "lambdar_disk.png")
    plt.close()    


#%%   
def cal_lambda_r(galaxy, npix=30, rscale=3.0, method=2, verbose=False):
    import numpy as np

    # Only within effective radius. (or whatever radius.)
    # already centered.
    rr = galaxy.reff * rscale
    
    ind_ok = np.where((abs(galaxy.star["x"]) <= rr) &
                      (abs(galaxy.star["y"]) <= rr) &
                      (abs(galaxy.star["z"]) <= rr) )[0]
# Some margin make mmap and vmap look better.
# But only inside Lambda_R for R < Reff is returned. 

    print("{} particles inside {}kpc".format(len(ind_ok), 
                rr))

    x = galaxy.star['x'][ind_ok]
    y = galaxy.star['y'][ind_ok]
    m = galaxy.star['m'][ind_ok]
    vz = galaxy.star['vz'][ind_ok]
    
    # NGP charge assignment         
    nx = npix
    ny = npix

    # fix center explicitely.
    # 0.5 * (min + max) != center
    xx = (x + rr) / rr * 0.5 * nx
    yy = (y + rr) / rr * 0.5 * ny
    
    # coordiantes of particles.
    ngx = np.clip(np.fix(xx), 0, nx-1)
    ngy = np.clip(np.fix(yy), 0, ny-1)
    
    indices = ngx + ngy * nx

    mmap = np.zeros(nx * ny, dtype=float)
    for i, ind in enumerate(indices):
        mmap[ind] += m[i]
    
    dx = (2 * rr) / npix
    mmap = mmap / (dx*dx)
    galaxy.mmap = mmap.reshape(nx, ny)
    
    vmap = np.zeros(nx * ny, dtype=float)
    # mass-weighted sigma         
    sigmap=np.zeros(nx * ny, dtype=float)
    for i in range(nx * ny):
        ind = np.where(indices == i)[0]
        if len(ind) > 0:
            sigmap[i] = galaxy.weighted_std(vz[ind], m[ind])
            vmap[i] = np.average(vz[ind], weights=m[ind])
        else:
            sigmap[i] = 0
            vmap[i] = 0

    points = np.zeros(0.5 * npix) # 1.5 Reff
    npoints = len(points)
    # distance in pixel unit.
    dist2d=np.zeros((nx, ny), dtype=float)
    for i in range(nx):
        for j in range(ny):
            dist2d[i][j]= np.sqrt((0.5 + i - nx/2)**2 + (0.5 + j - ny/2)**2)
            
    dist1d = dist2d.ravel()


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
                sig = galaxy.weighted_std(vz[ind], m[ind])
                b = sum(m[ind] * dd[ind] * np.sqrt(vv**2 + sig**2))
                points[i] = a/b
        
    lambda_arr = points
#    lambda_r = 0.5*(lambda_arr[0.5*len(lambda_arr) -1] 
#                            + lambda_arr[len(lambda_arr) -1])
    return dd, lambda_arr

#%%
def lambda_voronoi(galaxy, targetSN=15, plot=False, quiet=True,
                   npix=30, rscale=3.0, all_plot=True):
                       
    rr = galaxy.reff * rscale

    ind_ok = np.where((abs(galaxy.star["x"]) <= rr) &
                      (abs(galaxy.star["y"]) <= rr) &
                      (abs(galaxy.star["z"]) <= rr) )[0]
# Some margin make mmap and vmap look better.
# But only inside Lambda_R for R < Reff is returned. 

    print("{} particles inside {}kpc".format(len(ind_ok), 
                rr))

    x = galaxy.star['x'][ind_ok]
    y = galaxy.star['y'][ind_ok]
    m = galaxy.star['m'][ind_ok]
    vz = galaxy.star['vz'][ind_ok]
    
    # NGP charge assignment         
    nx = npix
    ny = npix

    # distance in pixel unit.
    dist2d=np.zeros((nx, ny), dtype=float)
    for i in range(nx):
        for j in range(ny):
            dist2d[i][j]= np.sqrt((0.5 + i - nx/2)**2 + (0.5 + j - ny/2)**2)
#    dist_part = np.sqrt(x**2 + y**2)
    # fix center explicitely.
    # 0.5 * (min + max) != center
    xx = (x + rr) / rr * 0.5 * nx
    yy = (y + rr) / rr * 0.5 * ny
    
    ngx = np.clip(np.fix(xx).astype(int), 0, nx-1)
    ngy = np.clip(np.fix(yy).astype(int), 0, ny-1)
    
    indices = ngx + ngy * nx # npart-length array             
   
    # total mass in each grid.
    # the lambdar calculation here is mass-weighted.
    mmap = np.zeros(nx * ny, dtype=float)
    for i, ind in enumerate(indices):
        mmap[ind] += m[i]
    
    dx = (2 * rr) / npix
    mmap = mmap / (dx*dx)
    
    vmap = np.zeros(nx * ny, dtype=float)
    # mass-weighted sigma         
    sigmap=np.zeros(nx * ny, dtype=float)
    for i in range(nx * ny):
        ind = np.where(indices == i)[0]
        if len(ind) > 0:
            sigmap[i] = galaxy.weighted_std(vz[ind], m[ind])
            vmap[i] = np.average(vz[ind], weights=m[ind])
        else:
            sigmap[i] = 0
            vmap[i] = 0    
    
    # In reality, mmap = count map
    count_map = np.zeros(nx * ny, dtype=float) 
    for i in range(nx * ny):
        ind = np.where(indices == i)[0]
        if len(ind) > 0:
            count_map[i] = len(ind)

    i_ok = np.unique(indices)
    count_map = count_map[i_ok]
    noise_map = np.sqrt(count_map)

    xpos_regular = np.zeros((nx, ny), dtype=float)
    ypos_regular = np.zeros((nx, ny), dtype=float)
    for i in range(nx):
        xpos_regular[i,:] = i
    for j in range(ny):
        ypos_regular[:,j] = j
    xpos_regular = xpos_regular.ravel()
    ypos_regular = ypos_regular.ravel()
    xpos_regular = xpos_regular[i_ok]
    ypos_regular = ypos_regular[i_ok]

    from Cappellari.voronoi.voronoi_2d_binning import voronoi_2d_binning
    """
    This function accepts only data on uniform grid...?
    """
    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = \
      voronoi_2d_binning(xpos_regular, ypos_regular, count_map,
                         noise_map, targetSN=targetSN,
                         plot=plot, quiet=True)
                                                 
    vmap_v = np.zeros(len(xNode)) # Not actually maps, but 1-D arrays.
    sigmap_v = np.zeros(len(xNode))
    lambdamap_v = np.zeros(len(xNode))                                                 

# Quantities in Voronio bin
    for ibin in np.arange(len(xNode)):
        ind = np.where(binNum == ibin)[0] # pixels in this Voronoi bin
        i_part = np.empty((0,0), dtype=int)
        for j in ind:
            i_part = np.append(i_part, np.where(indices == j)[0])
        # all particles belonging to one Voronoi cell
        vv = abs(np.average(vz[i_part], weights=m[i_part]))       
        sig = galaxy.weighted_std(vz[i_part], weights=m[i_part])
        vmap_v[ibin] = vv
        sigmap_v[ibin] = sig
    lambdamap_v = vmap_v / sigmap_v

#    _display_pixels(xpos_regular, ypos_regular, binNum, pixelSize=None)
# averaged values 
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
        vmap[ind] = vmap_v[ibin]
        sigmap[ind] = sigmap_v[ibin]

    if all_plot:
        plt.close()
        fig, axs = plt.subplots(2,2)
        fig.suptitle("ID: {}    z: {:2f}".format(str(galaxy.id).zfill(5),
                     galaxy.info.zred))

# Stellar particle density map
        ax = axs[0,0]
        ax.imshow(mmap.reshape(nx,ny),interpolation='nearest')
        
        axs[0,1].imshow(vmap.reshape(nx,ny),interpolation='nearest')
        axs[1,0].imshow(sigmap.reshape(nx,ny),interpolation='nearest')
        axs[1,1].plot(new_arr / new_cnt)
        axs[1,1].set_ylim(0,1)
        plt.show()

    # divide by mass.
    return new_arr / new_cnt


#%%------------------------------------------------------------------------
if __name__ == '__main__':

    wdir = '/home/hoseung/Work/data/05427/'
    nout=187
    # 691 : ~10000
    # 715 : 27518
    hh, info, star, dm, cell = prepare_data(base=wdir, nout=nout, idgal=1242)
    # Direct plot ---------------------------------------------------------
    plot_halo_region(hh.data, info, star, cell, npix=400, save_dir=wdir)
    
    galaxy = mk_gal(hh.data, info, save=False, rscale=0.3
                , verbose=False, save_dir=wdir,
                rscale_lambda=2.0, npix=400)
#%%
    import draw
    gas_map = draw.pp.pp_cell(cell, 100, info, verbose=False)
#%% 
    l_v = lambda_voronoi(galaxy, 15, npix=30, rscale=3.0,
                         plot=True, quiet=True, all_plot=True)

#%%
    rscale = 3.0
    npix=30
    plot=True
    quiet = True
    targetSN = 15
    rr = galaxy.reff * rscale

    ind_ok = np.where((abs(galaxy.star["x"]) <= rr) &
                      (abs(galaxy.star["y"]) <= rr) &
                      (abs(galaxy.star["z"]) <= rr) )[0]
# Some margin make mmap and vmap look better.
# But only inside Lambda_R for R < Reff is returned. 

    print("{} particles inside {}kpc".format(len(ind_ok), 
                rr))

    x = galaxy.star['x'][ind_ok]
    y = galaxy.star['y'][ind_ok]
    m = galaxy.star['m'][ind_ok]
    vz = galaxy.star['vz'][ind_ok]
    
    # NGP charge assignment         
    nx = npix
    ny = npix

    # distance in pixel unit.
    dist2d=np.zeros((nx, ny), dtype=float)
    for i in range(nx):
        for j in range(ny):
            dist2d[i][j]= np.sqrt((0.5 + i - nx/2)**2 + (0.5 + j - ny/2)**2)
    dist_part = np.sqrt(x**2 + y**2)
    # fix center explicitely.
    # 0.5 * (min + max) != center
    xx = (x + rr) / rr * 0.5 * nx
    yy = (y + rr) / rr * 0.5 * ny
    
    ngx = np.clip(np.fix(xx).astype(int), 0, nx-1)
    ngy = np.clip(np.fix(yy).astype(int), 0, ny-1)
    
    indices = ngx + ngy * nx # npart-length array             
   
    # total mass in each grid.
    # the lambdar calculation here is mass-weighted.
    mmap = np.zeros(nx * ny, dtype=float)
    for i, ind in enumerate(indices):
        mmap[ind] += m[i]
    
    dx = (2 * rr) / npix
    mmap = mmap / (dx*dx)
    
    vmap = np.zeros(nx * ny, dtype=float)
    # mass-weighted sigma         
    sigmap=np.zeros(nx * ny, dtype=float)
    for i in range(nx * ny):
        ind = np.where(indices == i)[0]
        if len(ind) > 0:
            sigmap[i] = galaxy.weighted_std(vz[ind], m[ind])
            vmap[i] = np.average(vz[ind], weights=m[ind])
        else:
            sigmap[i] = 0
            vmap[i] = 0    
    
    # In reality, mmap = count map
    count_map = np.zeros(nx * ny, dtype=float) 
    for i in range(nx * ny):
        ind = np.where(indices == i)[0]
        if len(ind) > 0:
            count_map[i] = len(ind)

    i_ok = np.unique(indices)
    count_map = count_map[i_ok]
    noise_map = np.sqrt(count_map)

    xpos_regular = np.zeros((nx, ny), dtype=float)
    ypos_regular = np.zeros((nx, ny), dtype=float)
    for i in range(nx):
        xpos_regular[i,:] = i
    for j in range(ny):
        ypos_regular[:,j] = j
    xpos_regular = xpos_regular.ravel()
    ypos_regular = ypos_regular.ravel()
    xpos_regular = xpos_regular[i_ok]
    ypos_regular = ypos_regular[i_ok]

#%%
    from Cappellari.voronoi.voronoi_2d_binning import voronoi_2d_binning
    """
    This function accepts only data on uniform grid...?
    """
    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = \
      voronoi_2d_binning(xpos_regular, ypos_regular, count_map,
                         noise_map, targetSN=targetSN,
                         plot=True, quiet=True)    

    vmap_v = np.zeros(len(xNode)) # Not actually maps, but 1-D arrays.
    sigmap_v = np.zeros(len(xNode))
    lambdamap_v = np.zeros(len(xNode))                                                 

# Quantities in Voronio bin
    for ibin in np.arange(len(xNode)):
        ind = np.where(binNum == ibin)[0] # pixels in this Voronoi bin
        i_part = np.empty((0,0), dtype=int)
        for j in ind:
            i_part = np.append(i_part, np.where(indices == j)[0])
        # all particles belonging to one Voronoi cell
        vv = abs(np.average(vz[i_part], weights=m[i_part]))       
        sig = galaxy.weighted_std(vz[i_part], weights=m[i_part])
        vmap_v[ibin] = vv
        sigmap_v[ibin] = sig
    lambdamap_v = vmap_v / sigmap_v
    
    for ibin in range(len(xNode)):
        ind = np.where(binNum == ibin)[0]
        vmap[ind] = vmap_v[ibin]
        sigmap[ind] = sigmap_v[ibin]

#%%    
    # need galactocentric distances of Voronoi cells
    plt.close()
    fig = plt.figure()
    axes = fig.add_subplot(211)    
#    axes.plot(dd[dd.argsort()], lambdamap_v[dd.argsort()])
    axes.plot(l_v)
    axes.set_ylim(0,1)                                                 

#%%               
#    dd, ll = cal_lamdar_v(galaxy, npix=30, rscale=3.0, method="VCT", targetSN=15)

#    plt.close()
#    fig = plt.figure()
#    axes = fig.add_subplot(212)
#    axes.scatter(dd[dd.argsort()], ll[dd.argsort()])
#    axes.set_ylim(0,2)
    
#%%
#    lambdar = cal_lamdar(gal, npix=30, rscale=3.0, method="VCT", targetSN=4)
#    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = \
#         cal_lamdar(gal, npix=30, rscale=3.0, method="VCT", targetSN=4)
    # Calculate lambda_r -------------------------------------------------
    # output data to anohter class / or container 