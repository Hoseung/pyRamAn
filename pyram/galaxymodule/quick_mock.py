import numpy as np
from .load_sed import *

class Simplemock():
    """
    Calculate flux through a filter for all stellar particle.

    Example
    >>> from general import defaults
    >>> import quick_mock
    >>> import galaxymodule
    >>> dfl = defaults.Default()
    >>> nout = 312#782 # 312
    >>> s = load.sim.Sim(nout=nout)
    >>> gcat = tree.halomodule.Halo(nout=nout, is_gal=True)
    >>> gg = galaxymodule.rd_GM.Gal(nout, catalog=gcat.data[21].copy(), info=s.info)
    >>> gg.debug=False
    >>> make_gal.mk_gal(gg)
    >>> from galaxymodule import quick_mock as qmc
    >>> MockSED = qmc.Simplemock(repo=dfl.dir_repo+'sed/')
    >>> gg.star.Flux_u= MockSED.get_flux(star=gg.star, filter_name='u')
    >>> gg.star.Flux_g= MockSED.get_flux(star=gg.star, filter_name='g')
    >>> gg.star.Flux_r= MockSED.get_flux(star=gg.star, filter_name='r')
    >>> gg.star.Flux_i= MockSED.get_flux(star=gg.star, filter_name='i')
    >>> gg.star.Flux_z= MockSED.get_flux(star=gg.star, filter_name='z')
    >>> Fluxs = [gg.star.Flux_u,
                 gg.star.Flux_g,
                 gg.star.Flux_r,
                 gg.star.Flux_i,
                 gg.star.Flux_z]
    >>> quick_mock.draw(Fluxs, gg.star["x"], gg.star["y"], suffix="face")
    >>> quick_mock.draw(Fluxs, gg.star["x"], gg.star["z"], suffix="edge")

    TODO

    t_univ

    interpolation outside range. (metallcity < 0.0004)

    """
    def __init__(self, repo="/home/hoseung/Work/pyclusterevol/repo/sed/",
                 filter_system="SDSS",
                 sed_model="bc03yy",
                 info=None,
                 load=True,
                 imf = "Salpeter"):
        self.filter_system = filter_system
        self.repo = repo
        self.sed_model = sed_model
        self.IMF = imf
        self.info = info
        if load is True:
            self.load_filters()
            self.load_SED_wavelength()
            self.load_SED_all()

    def load_filters(self):
        if self.filter_system == "SDSS":
            filter_lambda, filter_u, filter_g, filter_r, filter_i, filter_z = \
                np.genfromtxt(self.repo + "filter_sdss.dat",
                                skip_header=1, unpack=True)

            self.filters = {"lambda":filter_lambda,
                            "u":filter_u,
                            "g":filter_g,
                            "r":filter_r,
                            "i":filter_i,
                            "z":filter_z}

    def load_SED_wavelength(self):
        if self.sed_model == "bc03yy":
            self.sed_wavelength = np.genfromtxt(self.repo + "lambda.dat")

    def load_SED_all(self):
        """
        Full SEDs are just a few tens of MB.
        """
        if self.sed_model == "bc03yy":
            self.metal_points = np.array([0.0004, 0.001, 0.004, 0.01, 0.02, 0.04])
            # age points in tables.
            self.age_points = np.genfromtxt(self.repo+"ages_yybc.dat") # in Gry unit
            self.SEDs = np.zeros((6, 221, 1221))

            for i, metal in enumerate(self.metal_points):
                self.SEDs[i,:,:] = np.genfromtxt(self.repo +
                                            "bc03_yy_{:.4f}".format(metal)).reshape(221, 1221)
        elif self.sed_model == "bc":
            all_sed=[]
            for fn in glob("./BC03_v03/bc03/models/Padova2000/salpeter/*lr*.ised*.gz"):
                bc03 = Sed_block(fn)
                all_sed.append(bc03)
                print(bc03.isochrone)
                print(bc03.X, bc03.Y, bc03.Z)

            SED = Sed(all_sed)
            self.SEDs

        else:
            print("Sorry, Only bc03 is implemented.")

    def get_flux(self, star,
                 cell=None,
                 simple=False,
                 metal_lower_cut = True,
                 filter_name='r',
                 quick=False,
                 speed_check=False):
        """
        calculate SED of each particle.
        If cell is not None, attenuate flux accordingly.

        parameters
        ----------
        quick : False
            If True, reduce wavelength point of sed and filter significantly.
            This still gives reasonable value.

        speed_check : False
            If True, measure the time taken. Only for optimizing purpose.
        """
        if speed_check:
            from time import time
            t0 = time()

        Lum_sun = 3.826e33
        # BC2003 is in unit of L_sun Ang-1, where L_sun = Lum_sun.

        starmetal = star["metal"].copy() # Is the original array modified?
        starmetal[starmetal > 0.04] = 0.0399999 # a few stars have higher metallicity
        if metal_lower_cut:
            # No star with metallicity lower than the lowest table.
            starmetal[starmetal < min(self.metal_points)] = min(self.metal_points) * 1.0001

        locate_metal = np.digitize(starmetal, self.metal_points)-1 # GOOD
        relevant_metals = self.metal_points[:max(locate_metal)+2]
        nmetals = len(relevant_metals)

        # Star Age
        starage = star["age"] # the field name "time" should have been changed to "age"

        locate_age = np.digitize(starage, self.age_points)-1 # GOOD
        relevant_ages = self.age_points[:max(locate_age)+2]
        nages = len(relevant_ages)
        if speed_check: t1 = time() #

        ### Filter optimization. #################################################
        # Pick one
        this_filter = self.filters[filter_name]

        # band range
        i_filter_pos = this_filter > 0
        if quick:
            this_filter = this_filter[i_filter_pos][::10]
            filter_lambda_this_band = self.filters["lambda"][i_filter_pos][::10]
        else:
            this_filter = this_filter[i_filter_pos]
            filter_lambda_this_band = self.filters["lambda"][i_filter_pos]

        lambda_min_this_band = min(filter_lambda_this_band)
        lambda_max_this_band = max(filter_lambda_this_band)

        if quick:
            n_compress = 20
            sed_org = self.sed_wavelength
            sed_wavelength = self.sed_wavelength[::n_compress]
        else:
            sed_wavelength = self.sed_wavelength

        i_lambda_min = np.argmax(sed_wavelength > lambda_min_this_band) -1
        i_lambda_max = np.argmax(sed_wavelength > lambda_max_this_band)

        # Only a small part of SED is needed.
        # To compute d_lambda, one additional lambda point is desired.
        # Could be forward / backward / midpoint and so on.
        # let me take backward as fractional chnge in d_lambda is less in longer wavelength
        # Well.. actually I don't care..
        # d_lambda = wavelength[:-1] - wavelength[1:]
        wavelength = sed_wavelength[i_lambda_min:i_lambda_max+2] # why +2?
        n_wavelength = len(wavelength)-1#i_lambda_max - i_lambda_min + 1

        ##### Caclulate band flux #################
        # Load only necessary data
        # Load all once, keep under the class and copy a part of it when needed here.
        seds = np.zeros((nmetals, nages, n_wavelength)) # metal age lambda
        if self.sed_model == "bc03yy":
            for i, metal in enumerate(relevant_metals):
                if quick:
                    for j in range(seds.shape[1]):
                        seds[i,j,:] = np.interp(wavelength[:-1],
                                            sed_org,
                                            self.SEDs[i,j,:])
                else:
                    seds[i,:,:] = self.SEDs[i,:nages, i_lambda_min:i_lambda_max+1]

        if speed_check: t2 = time() # all set up

        # All are array calculations.
        # interpolation weight
        dl_m = (starmetal - relevant_metals[locate_metal] ) / \
                                     (relevant_metals[locate_metal+1] - relevant_metals[locate_metal])
        dr_m = (relevant_metals[locate_metal+1] - starmetal) / \
                                     (relevant_metals[locate_metal+1] - relevant_metals[locate_metal])
        dl_a = (starage - relevant_ages[locate_age] )   / \
                                     (relevant_ages[locate_age+1] - relevant_ages[locate_age])
        dr_a = (relevant_ages[locate_age+1] - starage ) / \
                                     (relevant_ages[locate_age+1] - relevant_ages[locate_age])

        if speed_check: t3 = time() # done first easy calculation

        # 2D linear interpolation
        # weight * SED.
        Flux =  np.multiply( (dr_m * dr_a), seds[locate_metal, locate_age,:].T).T +\
                np.multiply( (dl_m * dr_a), seds[locate_metal + 1, locate_age,:].T).T +\
                np.multiply( (dr_m * dl_a), seds[locate_metal, locate_age + 1, :].T).T +\
                np.multiply( (dl_m * dl_a), seds[locate_metal + 1, locate_age + 1,:].T).T

        if speed_check: t4 = time()
        # Convolve filter
        # Wavelengths at which filter function are defined are different from the SED wavelength points.
        # Interpolate filter function on SED points.
        filter_in_sed_wavelengths = np.interp(wavelength, filter_lambda_this_band, this_filter)
        Flux = np.multiply(filter_in_sed_wavelengths[:-1] * wavelength[-1], Flux)#\
        div = np.multiply(filter_in_sed_wavelengths[:-1], wavelength[-1])

        if speed_check: t5 = time()
        # Need to multiply stellar mass

        if cell is None or self.info is None or len(cell) == 0:
            return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"]
        else:
            #print(len(cell))
            colden = get_star_colden(star, cell) *self.info.unit_nH *self.info.unit_l / self.info.boxtokpc #* 1e3
            try:
                if colden == -1:
                    return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"]
            except:
                pass
            #lams = np.linspace(1e3,1e4,1e3) # in Angstrom
            waven = wavelength[:-1] * 1e-4 # in 1e-6 m

            Es = 0.44#*EBV
            # Hydrogen column number density
            colden = colden/5.8e21  # -1mag per 5.8e21 Hydrogen.
            tau = 0.4 * Es * np.outer(ext_curve_k(waven), colden)
            # colden = N_star array.
            # waven = N_wavelengths array
            # tau = N_wave X N_star
            # Flux = N_star X N_wave  ... (may be the other way around.)
            #print("tau", tau)
            F_ext= [ff*np.power(10,(-1*tt)) for ff, tt in zip(Flux, tau.T)]

            #print("before ext", np.sum(Flux, axis=1))
            #print("After ext", np.sum(F_ext, axis=1))
            if speed_check: 
                t6 = time()
                print("age metal digitize {:.3f}".format(t1-t0))
                print("setup done {:.3f}".format(t2-t0))
                print("interpolation coefficient {:.3f}".format(t3-t0))
                print("flux {:.3f}".format(t4-t0))
                print("filter, flux, div {:.3f}".format(t5-t0))
                print("After dust attanuation {:.3f}".format(t6-t0))
            return np.sum(F_ext, axis=1) / np.sum(div) * Lum_sun * star["m"]

##################################################################
def flux2mag(flux,
             gal,
             x1="x",
             x2="y",
             filter_name = 'r',
             gal_range=None,
             Lum_dist = 400,
             plate_scale = 0.24,
             npixmax = 1200):
    # Observation conndition
    # in Mpc.
    if gal_range is None:
        gal_range = [[-gal.meta.rgal,gal.meta.rgal]]*2

    npixx, npixy = get_npix(plate_scale, gal_range, Lum_dist, npixmax)

    # Calculate Unit.
    d_lum_10p = 3.0857e19 # lumminonsity distance of 10pc in cm
    speed_of_light = 3e18 # angstrom / sec
    kpc_to_cm = 3.0857e21
    ldcm = Lum_dist * kpc_to_cm
    inv_distance = 1/(4*np.pi * ldcm * ldcm)

    band=BandSDSS()
    # Additional factors to derive realistic flux values.

    print(npixx, npixy)
    Flux_map = np.histogram2d(gal.star[x1], gal.star[x2],
               weights=flux,
               bins=[npixx,npixy],
               range=gal_range)[0]

    Flux_map *= Lum_sun * 1e-2 * inv_distance
    return - 2.5 * np.log10(Flux_map) \
           - 5. * np.log10(getattr(band, filter_name)["pivot_lambda"]) \
           + 2.5 * np.log10(speed_of_light) -48.6


def get_npix(plate_scale, gal_range, Lum_dist, npixmax):
    FOVx = (gal_range[0][1] - gal_range[0][0]) / (Lum_dist*1e3) * 180. / np.pi * 3600. # in arcsec
    FOVy = (gal_range[1][1] - gal_range[1][0]) / (Lum_dist*1e3) * 180. / np.pi * 3600.
    npixx= int(np.ceil(FOVx/plate_scale))
    npixy= int(np.ceil(FOVy/plate_scale))

    npixx = min([npixmax, npixx])
    npixy = min([npixmax, npixy])

    return (npixx, npixy)

class BandSDSS():
    def __init__(self):
        self.u = dict(pivot_lambda = 3557.0, name="u")
        self.g = dict(pivot_lambda = 4702.0, name="g")
        self.r = dict(pivot_lambda = 6175.0, name="r")
        self.i = dict(pivot_lambda = 7491.0, name="i")
        self.z = dict(pivot_lambda = 8946.0, name="z")

def composite_rgb(x,y, weight_r, weight_g, weight_b,
                  npix=100,
                  multiply_r=1.3,
                  multiply_g=1.1,
                  multiply_b=1.0,
                  log_scale = True,
                  range=None):
    import numpy as np
    from PIL import Image

    rgbArray = np.zeros((npix,npix,3))
    rgbArray[..., 0] = np.histogram2d(x, y,
               weights=weight_r, bins=npix, range=range)[0] * multiply_r
    rgbArray[..., 1] = np.histogram2d(x,y,
               weights=weight_g, bins=npix, range=range)[0] * multiply_g
    rgbArray[..., 2] = np.histogram2d(x,y,
               weights=weight_b, bins=npix, range=range)[0] * multiply_b

    if log_scale:
        rgbArray = np.log10(rgbArray+1)
    rgbArray = rgbArray/rgbArray.max() * 255
    #return rgbArray

    img = Image.fromarray(rgbArray.astype('uint8'))

    return img.rotate(90)
    #img.save('myimg.jpeg')

def draw(gal,
         x1="x",
         x2="y",
         suffix="edge",
         npix=200,
         gal_range=None,
         R="Flux_g",
         G="Flux_r",
         B="Flux_i",
         cr=1.0, cg=1.0, cb=3.0):
    """
    Parameers
    ---------
    gal:

    suffix:

    npix:
        default = 200
    gal_range:
        image 2d span in kpc.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    channel_r = getattr(gal.star, R)
    channel_g = getattr(gal.star, G)
    channel_b = getattr(gal.star, B)

    if gal_range is None:
        gal_range = [[-gal.meta.rgal,gal.meta.rgal]]*2

    fig, axs = plt.subplots(2,3)
    fig.set_size_inches(8,6)

    titles=["u","g","r","i","z","composite"]
    for i, ax in enumerate(axs.ravel()[:5]):
        ax.hist2d(gal.star[x1], gal.star[x2],
               weights=np.log10(getattr(gal.star, "Flux_"+titles[i])+1),
               bins=npix,
               cmap=plt.cm.binary_r,
               norm=LogNorm(),
               range=gal_range)
        ax.set_aspect("equal")
        ax.set_title(titles[i])

    # Try scaling Flux_x arrays to make a better composite image
    comp_img = composite_rgb(gal.star[x1], gal.star[x2],
                             channel_r,
                             channel_g,
                             channel_b,
                             npix=npix,
                             multiply_r = cr,
                             multiply_g = cg,
                             multiply_b = cb,
                             log_scale=True,
                             range=gal_range)
    ax = axs[-1][-1]
    ax.imshow(comp_img)
    ax.set_title("composite")
    labels = [item.get_text() for item in ax.get_xticklabels()]
    empty_string_labels = ['']*len(labels)
    ax.set_aspect("equal")
    ax.set_xticklabels(empty_string_labels)
    ax.set_yticklabels(empty_string_labels)

    plt.savefig(str(gal.meta.id).zfill(5) + suffix + ".png", dpi=200)


def get_absolute_mag(flux, band=None, bandname="r"):
    """
    Don't forget to pass band.
    Initializing a new band class is expensive.

    Todo
    ----
    Figure out, why 1e-2??
    """
    if band is None:
        band = BandSDSS()

    speed_of_light = 3e18 # angstrom / sec
    d_lum_10p = 3.0857e19 # lumminonsity distance of 10pc in cm
    return - 2.5 * np.log10(np.sum(flux)/(4.*np.pi*d_lum_10p*d_lum_10p)*1e-2) \
            - 5. * np.log10(getattr(band, bandname)["pivot_lambda"]) \
           + 2.5 * np.log10(speed_of_light) - 48.6


def get_star_colden(star, cell):
    nstar = len(star)
    # no need to perform calculations on irrelevant cells
    # cut cells behind the further
    # This shrinks the cell length significantly.

    cell_dx_min = np.min(cell["dx"])

    ddx_h = cell["dx"] * 0.5

    sub_cell = cell[((cell["x"]+ddx_h) > star["x"].min()) * (cell["x"]-ddx_h < star["x"].max()) *\
                    ((cell["y"]+ddx_h) > star["y"].min()) * (cell["y"]-ddx_h < star["y"].max()) *\
                    ((cell["z"]+ddx_h) < star["z"].max())]

    #print("len subcell", len(sub_cell))
    if len(sub_cell) < 1:
        return -1

    ddx_h = sub_cell["dx"] * 0.5
    xl = sub_cell["x"] - ddx_h
    xr = sub_cell["x"] + ddx_h
    yl = sub_cell["y"] - ddx_h
    yr = sub_cell["y"] + ddx_h

    xrange = (xl.min(),xr.max())
    yrange = (yl.min(),yr.max())
    xspan = xrange[1]-xrange[0]
    yspan = yrange[1]-yrange[0]
    npixx = np.int(np.ceil(xspan/cell_dx_min))
    npixy = np.int(np.ceil(yspan/cell_dx_min))

    h = np.histogram2d(star["x"], star["y"],
                       bins=[npixx,npixy],
                       range=[xrange, yrange])

    dxmap = dymap = cell_dx_min

    # Sort cells
    i_relev = np.where(h[0].ravel() > 0)[0]

    # Sort stars
    ix=np.searchsorted(h[1], star["x"]) - 1
    iy=np.searchsorted(h[2], star["y"]) - 1

    ixy = iy + npixx*ix
    #print("histogram indicies and ixy are equivalent:", np.all(np.unique(ixy) == i_relev))

    i_star_sort = np.argsort(ixy)
    i_star = np.concatenate((np.searchsorted(ixy[i_star_sort], np.unique(ixy)),
                             [nstar]))

    # stars grouped in each lum map.
    sorted_star = star[i_star_sort]
    #sorted_flux = sorted_star["Flux_r"]#[i_star_sort]
    colden_star = np.zeros(nstar)

    #print('len star', len(i_star))
    #print("len cell", len(sub_cell))
    #print("minx cell", min(xl))
    #print("maxx cell", max(xr))
    #print("miny cell", min(yl))
    #print("maxy cell", max(yr))

    for i in range(len(i_star)-1):
        stars_here = sorted_star[i_star[i]:i_star[i+1]]
        jx, jy = ix[i_star_sort[i_star[i]]], iy[i_star_sort[i_star[i]]]
   #     print(jx, jy)
   #     print(h[1][jx], h[2][jy])
        cells_here = sub_cell[ (xr >= h[1][jx]) * (xl <= h[1][jx]) *\
                               (yr >= h[2][jy]) * (yl <= h[2][jy]) ]
   #     print(len(cells_here), "cells here")
        if len(cells_here) == 0:
            # at the edge
            colden_star[i_star[i]:i_star[i+1]] = 0
        else:
            # column density for each star.
            i_star_z = np.searchsorted(cells_here["z"]+cells_here["dx"], stars_here["z"])
            i_star_z[i_star_z >= len(cells_here)] = len(cells_here)-1
            colden_star[i_star[i]:i_star[i+1]] = np.cumsum(cells_here["var0"]*cells_here["dx"])[i_star_z]

    return colden_star


def ext_curve_k(lam, Rv=4.05):
    H_frac = 0.76

    lambda1 = 0.48613 # 1e-6 m
    lambda2 = 0.65628 # 1e-6 m

    # Why use only lambda 2
    inv_lambda1 = 1./lam[lam < lambda2]
    inv_lambda2 = 1./lam[lam > lambda2]

    k = np.zeros(len(lam))
    # lambda : 0.09 ~ 0.63 1e-6m
    k1 = 2.659 * (-2.156 + (1.509 * inv_lambda1)
                      - (0.198*inv_lambda1**2)
                  + (0.011*inv_lambda1**3)) + Rv
    # lambda : 0.63 ~ 5.08 1e-6m
    k2 = 2.659 * (-1.857 + (1.040*inv_lambda2)) + Rv

    k[lam < lambda2] = k1
    k[lam > lambda2] = k2
    return k
