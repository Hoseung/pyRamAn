# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 21:50:05 2015

Defines classes for plot.
2d density map, 3D particle distribuion and so on.
Along with the data, some meta data should accompany.


@author: hoseung
"""
class MapImg():
    """
    elementary wrapper of single 2D image.
    """
    def __init__(self, ptype=None, proj=None, npix=800, info=None, data=None):
        self.set_data(data)
        if proj is not None:
            self.set_proj(proj)
        if self.data is not None:
            self.npixx = self.data.shape[0]
            self.npixy = self.data.shape[1]
        if info is not None:
            self.pboxsize = info.pboxsize

    def set_region(self, region):
        self.region = region

    def set_npix(self, npix):
        self.npix = npix

    def set_ptype(self, ptype):
        """Why do you need ptyep???"""
        self.ptype = ptype

    def set_data(self, data=None):
        print("Updating data")
        self.data = data

    def set_proj(self, proj):
        self.proj = proj
        if self.proj == 'x':
            self.xr = "yr"
            self.yr = "zr"
            self.xlabel = "y position [Mpc]"
            self.ylabel = "z position [Mpc]"
        if self.proj == 'y':
            self.xr = "xr"
            self.yr = "zr"
            self.xlabel = "x position [Mpc]"
            self.ylabel = "z position [Mpc]"
        if self.proj == 'z':
            self.xr = "xr"
            self.yr = "yr"
            self.xlabel = "x position [Mpc]"
            self.ylabel = "y position [Mpc]"

        if self.xr == 'xr':
            self.xlabel = "x position [Mpc]"

    def show_data(self):
        """ Simple projection of the data
        """
        import matplotlib.pyplot as plt
        plt.imshow(self.data)
        plt.show()

    def pickle_data(self, fout=None):
        import pickle
        if fout is None:
            pass
        with open(fout, 'wb') as f:
            pickle.dump(self, f)

    def cal_field(self, part, npix=800):
        from draw import pp
        x = part.dm['px']
        y = part.dm['py']  # These are views. right?
        z = part.dm['py']
        m = part.dm['m'] * self.info.msun # part must be normalized already!
        self.set_data(data = pp.den2d(x, y, z, m, npix, self.info, cic=True, norm_integer=True))

    def plot_2d_den(self, ax=None, vmin=None, vmax=None, cname='brg', save=False,
                     show=True, dpi=100, nticks=3, zposition=True, log=True,
                     cbar_name=r'$Log\Sigma$ $(M_{\odot} kpc^{-2})$', **kwargs):
        """
        plot 2d image array.
        Works with both particle map and hydro map.
        But hydro map requires vmax and vmin to be explicitely given.
        
        -------
        To do: give reasonable vmin and vmax also for hydro map
        
        """
        import matplotlib.pyplot as plt
        import numpy as np
        import utils.prettyplot as ptt

        if vmax is None:
            vmax = self.data.max()

        if vmin is None:
            vmin = vmax * 1e-6

        if not show:
             plt.ioff()

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if log:
            from matplotlib.colors import LogNorm
            p = ax.imshow(self.data, cmap=plt.get_cmap(cname),
                          norm=LogNorm(vmin=vmin, vmax=vmax), **kwargs)
        else:
            p = ax.imshow(self.data, cmap=plt.get_cmap(cname), **kwargs)

        #print(self.pboxsize)         print(self.region['zc'])
        if self.region is not None:
            xr = np.asarray(self.region['xr'])
            yr = np.asarray(self.region['yr'])
            # region['xr'] and region['yr'] are first and second axes.
            # Not absolute axes.

            x_old, x_new = ptt.tickInPboxInt(xr * self.pboxsize, self.data.shape[0], nticks=nticks)
            y_old, y_new = ptt.tickInPboxInt(yr * self.pboxsize, self.data.shape[1], nticks=nticks)

#           axes.set_xlabel(self.xlabel) # = xr or yr or zr depending on self.proj.
#            axes.set_ylabel(self.ylabel)
            ax.set_xticks(x_old * self.data.shape[0])
            ax.set_xticklabels(x_new)
            ax.set_yticks(y_old * self.data.shape[1])
            ax.set_yticklabels(y_new)
            if zposition:
                annotation = 'z position[Mpc/h]: {:.2f}'.format(self.region["zc"] * self.pboxsize)
                ax.text(0.05, 1.02, annotation, transform = ax.transAxes, ha='left', fontsize=14)

        #legend
#        print("Min, Max value",vmin, vmax)
        cbar = plt.colorbar(p)
#        cbar.ax.set_yticklabels(['0','1','2','>3'])
        cbar.set_label(cbar_name, rotation=270, labelpad=10)

        if show:
            plt.show()
        if save:
            if save == True:
                print("No image name is given... default = output.png \n")
                save = 'output.png'
            plt.savefig(save, dpi=dpi)
            plt.close()


class MapSet(MapImg):
    """ Container of multiple images from the same simulation.
    May include iamges from mulitple nouts, well. I am not sure..


    Attributes
    ----------
    set_ptype :
        type of map ("dm", "star", "gas" and so on.. )
    info :
        general information of the simulation.
        nout, redshift, aexp, boxsize, etc.
    pboxsize, xr,yr,zr (or region)
    """
    def __init__(self, data=None,  info=None, region=None):
        pass

    def rename_attr(self, old, new):
        """
        Modifies name of attribute.
        name of old, and new attributes are required in strings.

        Parameters
        ----------

        old
            name of attribute
        new
            new name of the old attribute

        """
        setattr(self, new, getattr(self, old))
        delattr(self, old)

    def help(self):
        import pprint
        pprint.pprint(self.data)
        pprint.pprint(self.info)
        pprint.pprint(self.proj)

#   Other information may include the age of universe, particle type,