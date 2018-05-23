import numpy as np

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
