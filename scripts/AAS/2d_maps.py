
# coding: utf-8

# # 2d maps
# stellar/DM density maps, gas density/temperature/metallicity maps

# Related methods. 
# 
# pp.den2d, pp.part2den, draw.img,..
# 
# 
# def den2d(x, y, z, m, npix, region=None,
#                     ngp=False, cic=False, tsc=False,
#                     vmin=None, vmax=None, verbose=False,
#                     norm_integer=False):
#     """
#         accept view of x,y,z, arrays and return 2D projected density map.
# 
#         assignment funciton (ngp, cic, tsc) takes these arrays.
#         But no modification is required.
#     """
#     
#     
# 

# In[1]:
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import draw
import load 
import utils.sampling as smp
import tree.halomodule as hmo 
import numpy as np
import utils.util
# run time parameters
# 1. simulation parameter
wdir = './'

rscale = 2
npix = 1600
dpi = 250
extent = (0, npix, 0, npix)
#nouts = range(187, 188)
nouts = range(186,187)#[10, 30, 50, 80, 100, 150, 187]
#lmaxs = [19, 18, 18, 17, 17, 16, 16]
lmax = 17
hvar = "temp"
cbar=True
cmaps_heat = ["afmhot", "gist_heat", "CMRmap", "hot"][1:2]

plot_part = False
plot_hydro = True
plot_halo = False


hh = hmo.Halo(base=wdir, nout=26, halofinder='hm', is_gal=False)
#cluster = hh.data[np.argmax(hh.data['np'])]
cluster = hh.data[np.where(hh.data['np'] == max(hh.data["np"]))]
xc, yc, zc = cluster.x, cluster.y, cluster.z

for nout in nouts:
    info = load.info.Info(nout=nout, base=wdir, load=True)   
    if nout < 10:
        if nout < 7:
            tmax = 5.5
        else:
            tmax = None
        tmin = None
    else:         
        tmax = 7.8
        tmin = 4.2

    if nout < 26:
        do_region_search = False # in case of 05427
    else:
        do_region_search = True

    if do_region_search:
        # Find the central cluster.
        hh = hmo.Halo(base=wdir, nout=nout, halofinder='hm', is_gal=False, info=info, load=True)
        gg = hmo.Halo(base=wdir, nout=nout, halofinder='HM', is_gal=True,  info=info, load=True)
        i_center = np.where(hh.data['np'] == max(hh.data['np']))
        cluster = hh.data[i_center]
        xc, yc, zc, rr = cluster.x, cluster.y, cluster.z, cluster.rvir / info.aexp
        print(nout, xc,yc,zc,rr)
    #    continue
        
    #    xc,yc,zc,rr = 0.638454915, 0.76221734, 0.359773845, 0.00278157

        region = smp.set_region_multi(xc=xc, yc=yc, zc=zc, radius = rr * rscale)
    else:        
        rr = cluster.rvir / info.aexp
        region = smp.set_region_multi(xc=xc, yc=yc, zc=zc, radius = rr * rscale)


    s = load.sim.Sim(nout, base=wdir)
    s.set_ranges(region["ranges"])
    s.show_cpus()

    small_region = region

    if plot_part:
        ptypes=["star id pos mass", "dm id pos mass"]
        s.add_part(ptypes) 
        s.part.load(fortran=True)
        star = s.part.star
        dm = s.part.dm

        x = star['x']
        y = star['y']
        z = star['z']
        
        star_map = draw.pp.den2d(star['x'],star['y'],star['z'],star['m'] * info.msun / (info.pboxsize * info.pboxsize * 1e6)
                                 , npix, region=small_region
                                 , cic=True, norm_integer=False, proj='z', vmin=1000) # vmin in Msun / kpc2
        
        dm_map = draw.pp.den2d(dm['x'],dm['y'],dm['z'],dm['m'], npix, region=small_region, cic=True, norm_integer=False)

        ls = np.zeros((npix,npix))
        ii = star_map > 0
        ls[ii] = np.log10(star_map[ii]) # Stellar map HAS empty pixels.
        ls[~ii] = np.floor(ls.min())

        ii = dm_map > 0
        ld = np.zeros((npix,npix))
        ld[ii] = np.log10(dm_map[ii]) # Stellar map HAS empty pixels.
        ld[~ii] = np.floor(ld.min())


    if plot_hydro:
        s.add_hydro()
        s.hydro.amr2cell(lmax=lmax)
        cell = s.hydro.cell

    #small_centers = region['centers'].copy() # wihtout .copy(), reference is passed, and modifying small_ceters also modify region['centers']
    #small_centers[2] += 0.1 * region['radius'] # offset 
    #small_region = smp.set_region(centers=small_centers, radius=0.02 * rr[0])
        
        gas_map = draw.pp.pp_cell(cell, npix, s.info, region=small_region, verbose=True, hvar=hvar)   
    
        lg = np.log10(gas_map).copy() # gas_map does not have empty pixel.
        
    
    
    # In[176]:
    
    rgal = region['radius'] * s.info.pboxsize# * 1000
    
    if plot_part:
        ls[ls < 2] = np.nan
        fig, ax = plt.subplots()
        im2 = ax.imshow(ls, cmap='brg', interpolation='nearest', extent=extent, origin='lower')
        ax.set_xlabel("x [Mpc]")
        ax.set_xticks(np.linspace(0,npix,5))
        xticks = ["{:.1f}".format(x) for x in np.linspace(-rgal, rgal, num=5)]
        ax.set_xticklabels(xticks)
        
        ax.set_ylabel("y [Mpc]")
        ax.set_yticks(np.linspace(0,npix,5))
        yticks = ["{:.1f}".format(y) for y in np.linspace(-rgal, rgal, num=5)]
        ax.set_yticklabels(yticks)
        ax.set_title("z= {:.2f}".format(s.info.zred))
        if cbar:
            fig.colorbar(im2)
        
        plt.savefig(wdir + "2d_star_" + str(nout).zfill(3) + ".png", dpi=dpi)
        if plot_halo:
            draw.pp.pp_halo(gg, npix, region=small_region, rscale=1.0, ind=np.arange(len(gg.data)), axes=plt.gca(), name=True,
                    verbose=False, new_ax=False, linewidth=0.4, edgecolor='y')        
            plt.savefig(wdir + "2d_star_gal_id" + str(nout).zfill(3) + ".png", dpi=dpi)
        
            draw.pp.pp_halo(hh, npix, region=small_region, rscale=1.0, ind=i_center[0], axes=plt.gca(), name=False,
                    verbose=False, new_ax=False, linewidth=0.4, edgecolor='y')
        
            plt.savefig(wdir + "2d_star_halo_" + str(nout).zfill(3) + ".png", dpi=dpi)

        plt.close()
    
        ls = (ls - ls.min())/ls.ptp() * 256
        ls[ls < 100] = np.nan

    if plot_hydro:
        vmax = lg.max()
        vmin = lg.min() 
        #lg = (lg - lg.min())/lg.ptp() * 256 # normalize so that later better mixes with stellar / DM map. 

        for cmap in cmaps_heat:
            im3 = plt.imshow(np.transpose(lg), cmap=cmap, alpha=.8, interpolation='bilinear',
                             extent=extent, origin='lower', vmax=tmax, vmin=tmin)
            
            #im1 = plt.imshow(ls, cmap="Blues", interpolation='nearest', extent=extent,
            #                 origin='lower', alpha=.3)
            
            ax = plt.gca()
            ax.set_xlabel("x [Mpc]")
            ax.set_xticks(np.linspace(0,npix,5))
            xticks = ["{:.1f}".format(x) for x in np.linspace(-rgal, rgal, num=5)]
            ax.set_xticklabels(xticks)
            
            ax.set_ylabel("y [Mpc]")
            ax.set_yticks(np.linspace(0,npix,5))
            yticks = ["{:.1f}".format(y) for y in np.linspace(-rgal, rgal, num=5)]
            ax.set_yticklabels(yticks)
            ax.set_title("z= {:.2f}".format(s.info.zred))
            if cbar:
                plt.colorbar(im3, label="logT")
            
            plt.savefig(wdir + "2d_gas_" + str(nout).zfill(3) + hvar + cmap + ".png", dpi=dpi)
            
            if plot_halo:
                draw.pp.pp_halo(hh, npix, region=small_region, rscale=1.0, ind=i_center[0], axes=ax, name=False,
                        verbose=False, new_ax=False, linewidth=0.4)
                plt.savefig(wdir + "2d_gas_halo_" + str(nout).zfill(3) + hvar + cmap +".png", dpi=dpi)
            plt.close()



# By default(WHY????), (0,0) coordinate of an axes created by imshow starts from the top left. 
# origin='lower' to change the default coordinate origin. 

# In[88]:
    if plot_part:
        im1 = plt.imshow(ls, cmap="gray", interpolation='bilinear', extent=extent, origin='lower')
        im2 = plt.imshow(ld, cmap="copper", interpolation='bilinear', extent=extent, origin='lower')
        this_axes = plt.gca()
        hind = np.where(hh.data['mvir'] > 5e10)[0]
        if plot_halo:
            draw.pp.pp_halo(hh, npix, region=small_region, rscale=1.0, ind=hind, axes=this_axes, name=True,
                    verbose=False, new_ax=False)
        ax = plt.gca()    
        ax.set_title("z= {:.2f}".format(s.info.zred))
        plt.savefig(wdir + "2d_dm" + str(nout).zfill(3) +" .png", dpi=dpi)
        plt.close()
    
    
    # In[89]:
    if plot_part and plot_hydro:
        #plt.hold(True)
        im2 = plt.imshow(ld, cmap="copper", interpolation='nearest', extent=extent, origin='lower')
        im3 = plt.imshow(np.transpose(lg), cmap="CMRmap", alpha=.5, interpolation='bilinear', extent=extent, origin='lower')
        ax = plt.gca()    
        ax.set_title("z= {:.2f}".format(s.info.zred))
        plt.savefig(wdir + "2d_dm_gas" + str(nout).zfill(3) + hvar +".png", dpi=dpi)
        plt.close()
# One of two should be transposed.
# But which one?



# I don't understand.. 
# Anyways, small_region = smp.set_region(centers=small_centers, radius=0.05 * rr) results in small_region['centers'] to be a list of 3 arrays. Same goes for the other variables.
# But if you modify rr to rr[0], then the resulting region is good. 
# And, you need to pass the 'good' region to pp module in order to get the calculations right. 

# Note 1)
# When making an superimposed image, it is important to match the 'background color' of each figure.
# If low density gas has some color while empty stellar map is black, the resulting image has unwanted color that may indicate some other value of gas density, (or stellar density).
# 
# Note 2)
# region argument added to pp_cell.
# 
# Note 3)
# The above code does not run in parallel. 
# First suspect is "plt.hold(True)".
#   -> Do I really need plt.hold(True)?
