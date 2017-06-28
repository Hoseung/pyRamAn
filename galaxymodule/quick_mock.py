import numpy as np
class Simplemock():
    """
    Calculate flux through a filter for all stellar particle.

    Example
    >>> from general import defaults
    >>> import quick_mock
    >>> dfl = defaults.Default()
    >>> nout = 312#782 # 312
    >>> s = load.sim.Sim(nout=nout)
    >>> gcat = tree.halomodule.Halo(nout=nout, is_gal=True)
    >>> gg = load.rd_GM.Gal(nout, catalog=gcat.data[21].copy(), info=s.info)
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
                 sed_model="bc03",
                 load=True):
        self.filter_system = filter_system
        self.repo = repo
        self.sed_model = sed_model
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
        if self.sed_model == "bc03":
            self.sed_wavelength = np.genfromtxt(self.repo + "lambda.dat")

    def load_SED_all(self):
        """
        Full SEDs are just a few tens of MB.
        """
        if self.sed_model == "bc03":
            self.metal_points = np.array([0.0004, 0.001, 0.004, 0.01, 0.02, 0.04])
            # age points in tables.
            self.age_points = np.genfromtxt(self.repo+"ages_yybc.dat") # in Gry unit
            self.SEDs = np.zeros((6, 221, 1221))

            for i, metal in enumerate(self.metal_points):
                self.SEDs[i,:,:] = np.genfromtxt(self.repo +
                                            "bc03_yy_{:.4f}".format(metal)).reshape(221, 1221)
        else:
            print("Sorry, Only bc03 is implemented.")

    def get_flux(self, star,
                 extinction = False,
                 metal_lower_cut = True,
                 filter_name='r'):
        ### star data ########################################################
        # BC03 related.
        Lum_sun = 3.826e33
        # BC2003 is in unit of L_sun Ang-1, where L_sun = Lum_sun.


        starmetal = star["metal"] # Is the original array modified?
        if metal_lower_cut:
            # No star with metallicity lower than the lowest table.
            starmetal[starmetal < min(self.metal_points)] = min(self.metal_points) * 1.0001

        locate_metal = np.digitize(starmetal, self.metal_points)-1 # GOOD
        relevant_metals = self.metal_points[:max(locate_metal)+2]
        nmetals = len(relevant_metals)

        # Star Age
        starage = star["time"]

        locate_age = np.digitize(starage, self.age_points)-1 # GOOD
        relevant_ages = self.age_points[:max(locate_age)+2]
        nages = len(relevant_ages)

        ### Filter optimization. #################################################

        # Pick one
        this_filter = self.filters[filter_name]

        # band range
        i_filter_pos = this_filter > 0

        this_filter = this_filter[i_filter_pos]
        filter_lambda_this_band = self.filters["lambda"][i_filter_pos]

        lambda_min_this_band = min(filter_lambda_this_band)
        lambda_max_this_band = max(filter_lambda_this_band)

        i_lambda_min = np.argmax(self.sed_wavelength > lambda_min_this_band) -1
        i_lambda_max = np.argmax(self.sed_wavelength > lambda_max_this_band)
        #print(i_lambda_min, wavelength[i_lambda_min], lambda_min_this_band)
        #print(i_lambda_max, wavelength[i_lambda_max], lambda_max_this_band)

        # Only a small part of SED is needed.
        # To compute d_lambda, one additional lambda point is desired.
        # Could be forward / backward / midpoint and so on.
        # let me take backward as fractional chnge in d_lambda is less in longer wavelength
        # Well.. actually I don't care..
        # d_lambda = wavelength[:-1] - wavelength[1:]
        #
        wavelength = self.sed_wavelength[i_lambda_min:i_lambda_max+2]
        n_wavelength = i_lambda_max - i_lambda_min + 1

        ##### Caclulate band flux #################
        # Load only necessary data
        # Load all once, keep under the class and copy part of it when needed heere.
        seds = np.zeros((nmetals, nages, n_wavelength)) # metal age lambda
        if self.sed_model == "bc03":
            for i, metal in enumerate(relevant_metals):
                seds[i,:,:] = self.SEDs[i,:nages, i_lambda_min:i_lambda_max+1]

        # All are array-wise calculations.
        # interpolation weight
        dl_m = (starmetal - relevant_metals[locate_metal] ) / \
                                              (relevant_metals[locate_metal+1] - relevant_metals[locate_metal])
        dr_m = (relevant_metals[locate_metal+1] - starmetal) / \
                                              (relevant_metals[locate_metal+1] - relevant_metals[locate_metal])
        dl_a = (starage - relevant_ages[locate_age] )   / \
                                              (relevant_ages[locate_age+1] - relevant_ages[locate_age])
        dr_a = (relevant_ages[locate_age+1] - starage ) / \
                                              (relevant_ages[locate_age+1] - relevant_ages[locate_age])


        # 2D linear interpolation
        # weight * SED.
        Flux =  np.multiply( (dr_m * dr_a), seds[locate_metal, locate_age,:].T).T +\
                np.multiply( (dl_m * dr_a), seds[locate_metal + 1, locate_age,:].T).T +\
                np.multiply( (dr_m * dl_a), seds[locate_metal, locate_age + 1, :].T).T +\
                np.multiply( (dl_m * dl_a), seds[locate_metal + 1, locate_age + 1,:].T).T


        # Convolve filter
        # Wavelengths at which filter function are defined are different from the SED wavelength points.
        # Interpolate filter function on SED points.
        filter_in_sed_wavelengths = np.interp(wavelength, filter_lambda_this_band, this_filter)
        #d_lambda = filter_in_sed_wavelengths[1:] - filter_in_sed_wavelengths[:-1]
        #Flux = np.multiply(filter_in_sed_wavelengths[:-1] * \
        #                   (filter_in_sed_wavelengths[:-1] \
        #                    - filter_in_sed_wavelengths[1:]), Flux)
        Flux = np.multiply(filter_in_sed_wavelengths[:-1] * wavelength[-1], Flux)#\

        # Need to multiply stellar mass

        if not extinction:
            return np.sum(Flux, axis=1) * star["m"] * Lum_sun

        else:
            print("Extinction - Not yet implemented")
            return

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

#def magmap(x,y,flux):

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
