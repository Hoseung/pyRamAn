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
        #print(i_lambda_min, wavelength[i_lambda_min], lambda_min_this_band)
        i_lambda_max = np.argmax(self.sed_wavelength > lambda_max_this_band)
        #print(i_lambda_max, wavelength[i_lambda_max], lambda_max_this_band)

        # Only a small part of SED is needed.
        wavelength = self.sed_wavelength[i_lambda_min:i_lambda_max+1]
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
        Flux = np.multiply(filter_in_sed_wavelengths, Flux)

        if not extinction:
            return np.sum(Flux, axis=1)
        else:
            print("Extinction - Not yet implemented")
            return


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

def draw(Fluxs, x1, x2, suffix="edge", npix=200, cr=1.0, cg=1.0, cb=3.0):
    import numpy as np
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2,3)
    fig.set_size_inches(8,6)

    titles=["u","g","r","i","z","composite"]
    for i, ax in enumerate(axs.ravel()[:5]):
        flux = Fluxs[i]
        ax.hist2d(x1,x2,
               weights=np.log10(flux+1),
               bins=npix,
               cmap=plt.cm.binary_r,
               norm=LogNorm(),
               range=gal_range)
        ax.set_aspect("equal")
        ax.set_title(titles[i])

    # Try scaling Flux_x arrays to make a better composite image
    comp_img = composite_rgb(x1,x2,
                             Fluxs[3],
                             Fluxs[2],
                             Fluxs[1],
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

    plt.savefig(str(gg1.meta.id).zfill(5) + suffix + ".png", dpi=200)
