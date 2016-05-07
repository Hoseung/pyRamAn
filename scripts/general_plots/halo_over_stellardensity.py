# -*- coding: utf-8 -*-.
"""
Created on Thu May 28 22:10:36 2015

plot halos over stellar density map


@author: hoseung
"""
def plot_halo_region(region, info, star, cell,
                     npix=400,
                     fn_save='region_map.png',
                     show=True):
    import draw  
    extent = (0, npix, 0, npix)

    star_map = draw.pp.den2d(star['x'],star['y'],star['z'],star['m'], npix,
                             region=region, cic=True, norm_integer=False)
    if star_map is not False:
        ls = np.zeros((npix,npix))
        ii = star_map > 0
        ls[ii] = np.log10(star_map[ii]) # Stellar map HAS empty pixels.
        ls[star_map <= 0] = np.floor(ls.min())
        im1 = plt.imshow(ls, cmap="gray",
                         interpolation='nearest',
                         extent=extent,
                         origin='lower')
    
    # One of two should be transposed.
    # But which one?
    gas_map = draw.pp.pp_cell(cell, npix, info, region=region, verbose=False)
    im2 = plt.imshow(np.transpose(np.log10(gas_map)),
                     cmap="CMRmap",
                     alpha=.5,
                     interpolation='bilinear',
                     extent=extent,
                     origin='lower')

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
    
    plt.savefig(fn_save, dpi=144)
    if show:
        plt.show()
    else:
 #       plt.close()
        return

#%%

if __name__ == '__main__':
    #base = './'
    base = "./"
    frefine= 'refine_params.txt'
    fnml = 'cosmo_200.nml'
    #scale = 0.3
    ptype=["star pos mass"]
    hydro = True
    npix = 800
    dmo=False
    rscale = 1.0
    dpi = 200
    
    import load
    import utils.sampling as smp
    import draw
    import tree
    import numpy as np
    import matplotlib.pyplot as plt
    nout_ini = int(input("First snapshot: \n"))
    nout_end = int(input("Last snapshot: \n"))
    nouts = range(nout_ini, nout_end + 1)
    
    halofinder = "HM"
    #halofinder = "RS"
    #%%
    for nout in nouts:
        snout = str(nout).zfill(3)
        info = load.info.Info(nout=nout, base=base, load=True)

# Set region ------------------------------------
#        if halofinder == "HM":
        hall = tree.halomodule.Halo(nout=nout,
                                    base=base,
                                    halofinder=halofinder,
                                    info=info,
                                    load=True)

        i_center = np.where(hall.data['np'] == max(hall.data['np']))   
#            h_ind = smp.extract_halos_within(hall.data, i_center)
#            h = tree.halomodule.Halo()
#            h.derive_from(hall, ind=i_center[0])
        
        region = smp.set_region_multi(xc=hall.data['x'][i_center],
                                      yc=hall.data['y'][i_center],
                                      zc=hall.data['z'][i_center],
                                      radius = hall.data['rvir'][i_center]*rscale)
            
#%% Load data ------------------------------
        s = load.sim.Sim(nout, base, ranges=region["ranges"], setup=True)
        s.add_part(ptype, load=True, fortran=True)
        star = s.part.star
        s.add_hydro(lmax = 17, load=True)
        cell = s.hydro.cell
#%% draw plot and save it.
        fig = plt.figure()
        # is it OK to created a new figure object and not return it?
        axes = fig.add_subplot(111)  

        plot_halo_region(region, s.info, star, cell,
                         npix=npix, show=False)
 
#%%     
        i_hal_plot = np.where(hall.data['mvir'] > 1e11)[0]
   
        draw.pp.pp_halo(hall, npix, region=region,
                        rscale=1.0,
                        ind=i_hal_plot,
                        axes=None,
                        verbose=False,
                        name=True,
                        fontsize=8)
        plt.savefig(base + str(nout).zfill(3) + '.png', dpi=dpi)
        plt.close()