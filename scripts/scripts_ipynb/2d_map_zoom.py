
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
wdir = './29172/'

#nout = 307
#rscales = 0.5**(range(16))
#lmaxs = 20 - (8-np.repeat(range(8))
rscales = 0.5**(np.arange(10))
lmaxs = 19 - (8-np.arange(10))

npix = 512
dpi = 200
extent = (0, npix, 0, npix)
nouts = range(187, 188)
#nouts = range(186,187)#[10, 30, 50, 80, 100, 150, 187]
#lmaxs = [19, 18, 18, 17, 17, 16, 16]
hvars = ["rho", "temp", "metal"][1:2]
cbar=True
cmaps_heat = ["gist_heat", "viridis", "magma"]

plot_part = False
plot_hydro = True
plot_halo = False

for nout in nouts:
    info = load.info.Info(nout=nout, base=wdir, load=True)   
    tmax = None
    tmin = None

    hh = hmo.Halo(base=wdir, nout=nout, halofinder='hm', is_gal=False)
    #cluster = hh.data[np.argmax(hh.data['np'])]
    cluster = hh.data[np.where(hh.data['np'] == max(hh.data["np"]))]
    xc, yc, zc = cluster.x, cluster.y, cluster.z
    rr = cluster.rvir
    #rr = cluster.rvir / info.aexp

    for rscale_org, lmax in zip(rscales[8:], lmaxs[8:]):

        region = smp.set_region(xc=xc, yc=yc, zc=zc, radius = 0.5 * rscale_org)

        s = load.sim.Sim(nout, base=wdir)
        s.set_ranges(region["ranges"])
        s.show_cpus()

        s.add_hydro()
        lmax = min([19, lmax])
        s.hydro.amr2cell(lmax=lmax)
        cell = s.hydro.cell

        for ii in range(2):    
            # change rscale alone.
            # Do not read data again, 
            # no chnage in lmax.
            rscale = rscale_org * (1./np.sqrt(2))**(ii)
            print(" @@@@@@@ ")
            print("Lmax: ", lmax, "rscale: ", rscale)
            print(" @@@@@@@ \n\n")
            region = smp.set_region(xc=xc, yc=yc, zc=zc, radius = 0.5 * rscale)
            s.set_ranges(region["ranges"])

    #small_centers = region['centers'].copy() # wihtout .copy(), reference is passed, and modifying small_ceters also modify region['centers']
    #small_centers[2] += 0.1 * region['radius'] # offset 
            for hvar in hvars:
                gas_map = draw.pp.pp_cell(cell, npix, s.info, region=region,
                                          verbose=True, hvar=hvar, proj="y")
                if hvar == "none":
                    lg = gas_map.copy()
                else:
                    lg = np.log10(gas_map).copy() # gas_map does not have empty pixel.
                
                        
                rgal = region['radius'] * s.info.pboxsize# * 1000
                
                for cmap in cmaps_heat:
                    im3 = plt.imshow(np.transpose(lg), cmap=cmap, alpha=1.0, interpolation='bilinear',
                                     extent=extent, origin='lower', vmax=tmax, vmin=tmin)
                    
                    ax = plt.gca()
                    ax.set_xlabel("x [Mpc]")
                    ax.set_xticks(np.linspace(0,npix,5))
                    xticks = ["{:.2f}".format(x) for x in np.linspace(-rgal, rgal, num=5)]
                    ax.set_xticklabels(xticks)
                    
                    ax.set_ylabel("y [Mpc]")
                    ax.set_yticks(np.linspace(0,npix,5))
                    yticks = ["{:.2f}".format(y) for y in np.linspace(-rgal, rgal, num=5)]
                    ax.set_yticklabels(yticks)
                    ax.set_title("z= {:.2f}".format(s.info.zred))
                    if cbar:
                        if hvar == "temp":
                            plt.colorbar(im3, label="logT")
                        elif hvar == "rho":
                            plt.colorbar(im3, label="density")
                        elif hvar == "metal":
                            plt.colorbar(im3, label="metallicity")
                    
                    plt.savefig(wdir+"2d_gas_"+str(nout).zfill(3)+"_{:.4f}_".format(rscale)+hvar+"_"+cmap+".png", dpi=dpi)
                    
                    plt.close()

