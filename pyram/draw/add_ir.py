import numpy as np
import math
import scipy.integrate
#import measure_run_time as mrt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable

import load
import utils
import tree
import utils.sampling as smp
#from galaxymodule import quick_mock as qmc

import sys
import os

sys.path.append("/home/jangjk816/pycluesterevol/galaxymodule")
sys.path.append("/home/jangjk816/Project/Mock/week3")
sys.path
#from galaxymodule import quick_mock_original as qmc_ori
#import try_quick_mock as tqmc

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits


######### constants #####################3

kpc_in_cm = 3.0857e21
msun_in_g = 1.989e33

G = 6.674e-8  # gravitational constant ; [cm^3 g^-1 s^-2]

#cell_ddx_grid_size = [2.38418579e-07,4.76837158e-07,9.53674316e-07,1.90734863e-06,3.81469727e-06,7.62939453e-06]

mingrid = 2.38418579e-07
cell_ddx_grid_size = [(2**0)*mingrid,(2**1)*mingrid,(2**2)*mingrid,(2**3)*mingrid,
                      (2**4)*mingrid,(2**5)*mingrid,(2**6)*mingrid,(2**7)*mingrid]


######################################################



def cos(theta):
    # return np.cos( theta * np.pi / 180 )
    return np.cos(theta)


def sin(theta):
    # return np.sin( theta * np.pi / 180 )
    return np.sin(theta)


#### Rotation ####
## you don't have to seperate each planes. Just for convenience ##

def rotation_xy(theta, xarray, yarray):
    x = cos(theta) * xarray - sin(theta) * yarray
    y = sin(theta) * xarray + cos(theta) * yarray
    return x, y


def rotation_xz(theta, xarray, zarray):
    x = cos(theta) * xarray - sin(theta) * zarray
    z = sin(theta) * xarray + cos(theta) * zarray
    return x, z


def rotation_yz(theta, yarray, zarray):
    y = cos(theta) * yarray - sin(theta) * zarray
    z = sin(theta) * yarray + cos(theta) * zarray
    return y, z


#### Axis Transformation ####


def cartesian_to_cylinder(xarray, yarray, zarray, vxarray, vyarray, vzarray):
    r = (xarray ** 2 + yarray ** 2) ** 0.5
    phi = np.arctan2(yarray, xarray)
    z = zarray

    v_r = np.cos(phi) * vxarray + np.sin(phi) * vyarray
    v_phi = (-np.sin(phi) * vxarray + np.cos(phi) * vyarray)
    v_z = vzarray

    return r, phi, z, v_r, v_phi, v_z


def plane_fit(x, a, b, c):
    return a * x[0] + b * x[1] + c


def auto_rotation_np(param, xarray, yarray, zarray,thetype=0):
    a, b, c = param[0], param[1], param[2]
    theta_xy = -1 * np.arctan2(a, b)
    yarray_mid, xarray_new = rotation_xy(theta_xy, yarray, xarray)
    y_abs = (a ** 2 + b ** 2) ** 0.5
    theta_yz = -1 * np.arctan2(y_abs, c)
    zarray_new, yarray_new = rotation_yz(theta_yz, zarray, yarray_mid)

    print(theta_xy * 180 / np.pi, theta_yz * 180 / np.pi)
    if thetype == 1:
        return xarray_new, yarray_new, zarray_new, theta_xy * 180 / np.pi, theta_yz * 180 / np.pi
    else:
        return xarray_new, yarray_new, zarray_new

#### Draw ####

def gm2code(arr, info):
    return (arr / info.pboxsize + 0.5)



###############################################3

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

    def __init__(self, repo="/home/jangjk816/pyclusterevol/repo/sed/",
                 filter_system="SDSS",
                 sed_model="bc03",
                 info=None,
                 load=True):
        self.filter_system = filter_system
        self.repo = repo
        self.sed_model = sed_model
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

            self.filters = {"lambda": filter_lambda,
                            "u": filter_u,
                            "g": filter_g,
                            "r": filter_r,
                            "i": filter_i,
                            "z": filter_z}

        elif self.filter_system == "BVR":
            filter_lambda, filter_B, filter_V, filter_R = \
                    np.genfromtxt(self.repo + "filter_BVR.dat",
                                    skip_header=1, unpack=True)

            self.filters = {"lambda": filter_lambda,
                                    "B": filter_B,
                                    "V": filter_V,
                                    "R": filter_R}

        elif self.filter_system == "FUV":
            filter_lambda, filter_FUV = \
                    np.genfromtxt(self.repo + "filter_fuv.dat",
                                    unpack=True)

            self.filters = {"lambda": filter_lambda,
                                    "FUV": filter_FUV}

        elif self.filter_system == "JHK":
            filter_lambda, filter_J, filter_H, filter_K = \
                    np.genfromtxt(self.repo + "filter_JHK.dat",
                                    skip_header=1, unpack=True)

            self.filters = {"lambda": filter_lambda,
                                    "J": filter_J,
                                    "H": filter_H,
                                    "K": filter_K}



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
            self.age_points = np.genfromtxt(self.repo + "ages_yybc.dat")  # in Gry unit
            self.SEDs = np.zeros((6, 221, 1221))

            for i, metal in enumerate(self.metal_points):
                self.SEDs[i, :, :] = np.genfromtxt(self.repo +
                                                   "bc03_yy_{:.4f}".format(metal)).reshape(221, 1221)
        else:
            print("Sorry, Only bc03 is implemented.")

    def get_flux(self, star, num, colden=[],
                 extinction=False,
                 cell=None,
                 simple=False,
                 metal_lower_cut=True,
                 metal_higher_cut = True,
                 filter_name='r',
                 quick=False,
                 speed_check=False,
                 info = None):
        """
        calculate SED of each particle.
        If cell is not None, attenuate flux accordingly.

        parameters:
        ----------
        quick : False
            If True, reduce wavelength point of sed and filter significantly.
            This still gives reasonable value.

        speed_check : False
            If True, measure time taken. Only for optimizing purpose.
        """
        if speed_check:
            from time import time

        from time import time

        t0 = time()
        Lum_sun = 3.826e33
        # BC2003 is in unit of L_sun Ang-1, where L_sun = Lum_sun.

        starmetal = star["metal"].copy()  # Is the original array modified?
        starmetal[starmetal > 0.04] = 0.0399999  # a few stars have higher metallicity
        if metal_lower_cut:
            # No star with metallicity lower than the lowest table.
            starmetal[starmetal < min(self.metal_points)] = min(self.metal_points) * 1.0001
        if metal_higher_cut:
            # No star with metallicity lower than the lowest table.
            starmetal[starmetal > max(self.metal_points)] = max(self.metal_points) * 0.9999

        locate_metal = np.digitize(starmetal, self.metal_points) - 1  # GOOD
        relevant_metals = self.metal_points[:max(locate_metal) + 2]
        nmetals = len(relevant_metals)
        # print("Metal ranges", self.metal_points, starmetal.min(), starmetal.max())
        # print("nmetals", nmetals)

        # Star Age
        starage = -star["time"]

        tc = utils.cosmology.Timeconvert(info=info)

        starage = tc.time2gyr(star['time'],z_now = info.zred)

        locate_age = np.digitize(starage, self.age_points) - 1  # GOOD
        relevant_ages = self.age_points[:max(locate_age) + 2]
        nages = len(relevant_ages)
        t1 = time()  #

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

        i_lambda_min = np.argmax(sed_wavelength > lambda_min_this_band) - 1
        i_lambda_max = np.argmax(sed_wavelength > lambda_max_this_band)

        # Only a small part of SED is needed.
        # To compute d_lambda, one additional lambda point is desired.
        # Could be forward / backward / midpoint and so on.
        # let me take backward as fractional chnge in d_lambda is less in longer wavelength
        # Well.. actually I don't care..
        # d_lambda = wavelength[:-1] - wavelength[1:]
        wavelength = sed_wavelength[i_lambda_min:i_lambda_max + 2]  # why +2?
        n_wavelength = len(wavelength) - 1  # i_lambda_max - i_lambda_min + 1

        ##### Caclulate band flux #################
        # Load only necessary data
        # Load all once, keep under the class and copy part of it when needed heere.
        seds = np.zeros((nmetals, nages, n_wavelength))  # metal age lambda
        if self.sed_model == "bc03":
            for i, metal in enumerate(relevant_metals):
                if quick:
                    for j in range(seds.shape[1]):
                        seds[i, j, :] = np.interp(wavelength[:-1],
                                                  sed_org,
                                                  self.SEDs[i, j, :])
                else:
                    seds[i, :, :] = self.SEDs[i, :nages, i_lambda_min:i_lambda_max + 1]

        t2 = time()  # all set up

        # All are array calculations.
        # interpolation weight
        dl_m = (starmetal - relevant_metals[locate_metal]) / \
               (relevant_metals[locate_metal + 1] - relevant_metals[locate_metal])
        dr_m = (relevant_metals[locate_metal + 1] - starmetal) / \
               (relevant_metals[locate_metal + 1] - relevant_metals[locate_metal])
        dl_a = (starage - relevant_ages[locate_age]) / \
               (relevant_ages[locate_age + 1] - relevant_ages[locate_age])
        dr_a = (relevant_ages[locate_age + 1] - starage) / \
               (relevant_ages[locate_age + 1] - relevant_ages[locate_age])

        t3 = time()  # done first easy calculation

        # 2D linear interpolation
        # weight * SED.
        Flux = np.multiply((dr_m * dr_a), seds[locate_metal, locate_age, :].T).T + \
               np.multiply((dl_m * dr_a), seds[locate_metal + 1, locate_age, :].T).T + \
               np.multiply((dr_m * dl_a), seds[locate_metal, locate_age + 1, :].T).T + \
               np.multiply((dl_m * dl_a), seds[locate_metal + 1, locate_age + 1, :].T).T

        # print(np.sum(Flux, axis=1))

        t4 = time()
        # Convolve filter
        # Wavelengths at which filter function are defined are different from the SED wavelength points.
        # Interpolate filter function on SED points.
        filter_in_sed_wavelengths = np.interp(wavelength, filter_lambda_this_band, this_filter)
        Flux = np.multiply(filter_in_sed_wavelengths[:-1] * wavelength[-1], Flux)  # \
        div = np.multiply(filter_in_sed_wavelengths[:-1], wavelength[-1])
        print(wavelength)
        print(filter_lambda_this_band)
        print(this_filter)
        print(filter_in_sed_wavelengths[:-1])
        print(wavelength[-1])
        print(div)
        print(Flux,axis=1)


        t5 = time()
        # Need to multiply stellar mass
        # print(cell)
        # print(len(cell))
        # print(self.info)
        if extinction == True:
            if cell is None or self.info is None or len(cell) == 0:
                return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"] , []
            elif len(colden) ==0:
                return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"], []
            else :
                #print(Flux)

                # print(len(cell))
                #colden = get_star_colden(star, cell,
                #                         num) * self.info.unit_nH * self.info.unit_l / self.info.boxtokpc  # * 1e3
                if max(star['m']) < 1:
                    star["m"] *= 1e11

                print("Before extinction = %s" % np.sum((np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star['m'])))
                try:
                    if colden == -1:
                        return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"], []
                except:
                    pass
                # lams = np.linspace(1e3,1e4,1e3) # in Angstrom
                waven = wavelength[:-1] * 1e-4  # in 1e-6 m

                Es = 0.44  # *EBV
                # Hydrogen column number density
                # colden = colden/5.8e21  # -1mag per 5.8e21 Hydrogen.

                #print("colden is %s" % colden)
                #colden = colden * self.info.boxtokpc
                colden = colden / 5.8e21  # -1mag per 5.8e21 Hydrogen.
                tau = 0.4 * Es * np.outer(ext_curve_mw(waven), colden)
                #tau = 0.4 * Es * np.outer(ext_curve_k(waven), colden)
                # colden = N_star array.
                # waven = N_wavelengths array
                # tau = N_wave X N_star
                # Flux = N_star X N_wave  ... (may be the other way around.)
                # print("tau", tau)
                # print(np.sum(Flux,axis=1))

                F_ext = [ff * np.power(10, (-1 * tt)) for ff, tt in zip(Flux, tau.T)]

                # print(np.sum(F_ext, axis = 1))

                # print("before ext", np.sum(Flux, axis=1))
                # print("After ext", np.sum(F_ext, axis=1))
                t6 = time()
                # print("age metal digitize {:.3f}".format(t1-t0))
                # print("setup done {:.3f}".format(t2-t0))
                # print("interpolation coefficient {:.3f}".format(t3-t0))
                # print("flux {:.3f}".format(t4-t0))
                # print("filter, flux, div {:.3f}".format(t5-t0))
                # print("After dust attanuation {:.3f}".format(t6-t0))

                import matplotlib.pyplot as plt
                #print(len(F_ext[0]))
                forsed = []
                i = 0
                while i < len(tau):
                    forsed.append(tau[i][0])
                    i = i + 1
                # print(waven, tau)
                # if filter_name == 'u':
                # plt.plot(waven, Flux[0], c = 'black')
                # plt.plot(waven, F_ext[0], c = 'red')
                # plt.plot(waven, np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"], c = 'black')
                # plt.scatter(np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"],\

                # plt.legend()
                print("After Extinction = %s" % np.sum((np.sum(F_ext, axis=1) / np.sum(div) * Lum_sun * star['m'])))
                print('')


                return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"] , np.sum(F_ext, axis=1) / np.sum(div) * Lum_sun * star["m"]

        else :
            if cell is None or self.info is None or len(cell) == 0:
                return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"]
            elif len(colden) ==0:
                return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"]
            else :
                #print(Flux)

                # print(len(cell))
                #colden = get_star_colden(star, cell,
                #                         num) * self.info.unit_nH * self.info.unit_l / self.info.boxtokpc  # * 1e3
                if max(star['m']) < 1:
                    star["m"] *= 1e11

                print("Before extinction = %s" % np.sum((np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star['m'])))
                try:
                    if colden == -1:
                        return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"]
                except:
                    pass
                # lams = np.linspace(1e3,1e4,1e3) # in Angstrom
                waven = wavelength[:-1] * 1e-4  # in 1e-6 m

                Es = 0.44  # *EBV
                # Hydrogen column number density
                # colden = colden/5.8e21  # -1mag per 5.8e21 Hydrogen.

                #print("colden is %s" % colden)
                #colden = colden * self.info.boxtokpc
                colden = colden / 5.8e21  # -1mag per 5.8e21 Hydrogen.
                tau = 0.4 * Es * np.outer(ext_curve_k(waven), colden)
                # colden = N_star array.
                # waven = N_wavelengths array
                # tau = N_wave X N_star
                # Flux = N_star X N_wave  ... (may be the other way around.)
                # print("tau", tau)
                # print(np.sum(Flux,axis=1))

                F_ext = [ff * np.power(10, (-1 * tt)) for ff, tt in zip(Flux, tau.T)]

                # print(np.sum(F_ext, axis = 1))

                # print("before ext", np.sum(Flux, axis=1))
                # print("After ext", np.sum(F_ext, axis=1))
                t6 = time()
                # print("age metal digitize {:.3f}".format(t1-t0))
                # print("setup done {:.3f}".format(t2-t0))
                # print("interpolation coefficient {:.3f}".format(t3-t0))
                # print("flux {:.3f}".format(t4-t0))
                # print("filter, flux, div {:.3f}".format(t5-t0))
                # print("After dust attanuation {:.3f}".format(t6-t0))

                import matplotlib.pyplot as plt
                #print(len(F_ext[0]))
                forsed = []
                i = 0
                while i < len(tau):
                    forsed.append(tau[i][0])
                    i = i + 1
                # print(waven, tau)
                # if filter_name == 'u':
                # plt.plot(waven, Flux[0], c = 'black')
                # plt.plot(waven, F_ext[0], c = 'red')
                # plt.plot(waven, np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"], c = 'black')
                # plt.scatter(np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"],\

                # plt.legend()
                print("After Extinction = %s" % np.sum((np.sum(F_ext, axis=1) / np.sum(div) * Lum_sun * star['m'])))
                print('')
                return np.sum(F_ext, axis=1) / np.sum(div) * Lum_sun * star["m"]




    def get_SED(self, star, num, colden=[],
                 extinction=False,
                 cell=None,
                 simple=False,
                 metal_lower_cut=True,
                 metal_higher_cut = True,
                 filter_name='ext',
                 quick=False,
                 speed_check=False,
                 info = None):
        """
        calculate SED of each particle.
        If cell is not None, attenuate flux accordingly.

        parameters:
        ----------
        quick : False
            If True, reduce wavelength point of sed and filter significantly.
            This still gives reasonable value.

        speed_check : False
            If True, measure time taken. Only for optimizing purpose.
        """
        if speed_check:
            from time import time

        from time import time

        t0 = time()
        Lum_sun = 3.826e33
        # BC2003 is in unit of L_sun Ang-1, where L_sun = Lum_sun.

        starmetal = star["metal"].copy()  # Is the original array modified?
        starmetal[starmetal > 0.04] = 0.0399999  # a few stars have higher metallicity
        if metal_lower_cut:
            # No star with metallicity lower than the lowest table.
            starmetal[starmetal < min(self.metal_points)] = min(self.metal_points) * 1.0001
        if metal_higher_cut:
            # No star with metallicity lower than the lowest table.
            starmetal[starmetal > max(self.metal_points)] = max(self.metal_points) * 0.9999

        locate_metal = np.digitize(starmetal, self.metal_points) - 1  # GOOD
        relevant_metals = self.metal_points[:max(locate_metal) + 2]
        nmetals = len(relevant_metals)
        # print("Metal ranges", self.metal_points, starmetal.min(), starmetal.max())
        # print("nmetals", nmetals)

        # Star Age
        starage = -star["time"]

        tc = utils.cosmology.Timeconvert(info=info)

        starage = tc.time2gyr(star['time'],z_now = info.zred)

        locate_age = np.digitize(starage, self.age_points) - 1  # GOOD
        relevant_ages = self.age_points[:max(locate_age) + 2]
        nages = len(relevant_ages)
        t1 = time()  #

        ### Filter optimization. #################################################
        # Pick one
        filter_lambda, filter_ext = \
            np.genfromtxt(self.repo + "extinction_curve.dat",
                          skip_header=0, unpack=True)

        self.filters = {"lambda": filter_lambda,
                        "ext": filter_ext}


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

        i_lambda_min = np.argmax(sed_wavelength > lambda_min_this_band) - 1
        i_lambda_max = np.argmax(sed_wavelength > lambda_max_this_band)

        # Only a small part of SED is needed.
        # To compute d_lambda, one additional lambda point is desired.
        # Could be forward / backward / midpoint and so on.
        # let me take backward as fractional chnge in d_lambda is less in longer wavelength
        # Well.. actually I don't care..
        # d_lambda = wavelength[:-1] - wavelength[1:]
        wavelength = sed_wavelength[i_lambda_min:i_lambda_max + 2]  # why +2?
        n_wavelength = len(wavelength) - 1  # i_lambda_max - i_lambda_min + 1

        ##### Caclulate band flux #################
        # Load only necessary data
        # Load all once, keep under the class and copy part of it when needed heere.
        seds = np.zeros((nmetals, nages, n_wavelength))  # metal age lambda
        if self.sed_model == "bc03":
            for i, metal in enumerate(relevant_metals):
                if quick:
                    for j in range(seds.shape[1]):
                        seds[i, j, :] = np.interp(wavelength[:-1],
                                                  sed_org,
                                                  self.SEDs[i, j, :])
                else:
                    seds[i, :, :] = self.SEDs[i, :nages, i_lambda_min:i_lambda_max + 1]

        t2 = time()  # all set up

        # All are array calculations.
        # interpolation weight
        dl_m = (starmetal - relevant_metals[locate_metal]) / \
               (relevant_metals[locate_metal + 1] - relevant_metals[locate_metal])
        dr_m = (relevant_metals[locate_metal + 1] - starmetal) / \
               (relevant_metals[locate_metal + 1] - relevant_metals[locate_metal])
        dl_a = (starage - relevant_ages[locate_age]) / \
               (relevant_ages[locate_age + 1] - relevant_ages[locate_age])
        dr_a = (relevant_ages[locate_age + 1] - starage) / \
               (relevant_ages[locate_age + 1] - relevant_ages[locate_age])

        t3 = time()  # done first easy calculation

        # 2D linear interpolation
        # weight * SED.
        Flux = np.multiply((dr_m * dr_a), seds[locate_metal, locate_age, :].T).T + \
               np.multiply((dl_m * dr_a), seds[locate_metal + 1, locate_age, :].T).T + \
               np.multiply((dr_m * dl_a), seds[locate_metal, locate_age + 1, :].T).T + \
               np.multiply((dl_m * dl_a), seds[locate_metal + 1, locate_age + 1, :].T).T

        # print(np.sum(Flux, axis=1))

        t4 = time()
        # Convolve filter
        # Wavelengths at which filter function are defined are different from the SED wavelength points.
        # Interpolate filter function on SED points.
        filter_in_sed_wavelengths = np.interp(wavelength, filter_lambda_this_band, this_filter)
        Flux = np.multiply(filter_in_sed_wavelengths[:-1] * wavelength[-1], Flux)  # \
        div = np.multiply(filter_in_sed_wavelengths[:-1], wavelength[-1])

        t5 = time()
        # Need to multiply stellar mass
        # print(cell)
        # print(len(cell))
        # print(self.info)

        #print(Flux)

        # print(len(cell))
        #colden = get_star_colden(star, cell,
        #                         num) * self.info.unit_nH * self.info.unit_l / self.info.boxtokpc  # * 1e3
        if max(star['m']) < 1:
            star["m"] *= 1e11

        
        Flex = np.zeros(len(Flux[0]))
        l = 0
        while l < len(Flux):
            m = 0
            
            Flex += Flux[l] * star['m'][l]
            l += 1
            print(l/len(Flux)*100)

        print(Flex)


        print("Before extinction = %s" % np.sum((np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star['m'])))
        try:
            if colden == -1:
                return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"], []
        except:
            pass
        # lams = np.linspace(1e3,1e4,1e3) # in Angstrom
        waven = wavelength[:-1] * 1e-4  # in 1e-6 m

        Es = 0.44  # *EBV
        # Hydrogen column number density
        # colden = colden/5.8e21  # -1mag per 5.8e21 Hydrogen.

        #print("colden is %s" % colden)
        #colden = colden * self.info.boxtokpc
        colden = colden / 5.8e21  # -1mag per 5.8e21 Hydrogen.
        tau = 0.4 * Es * np.outer(ext_curve_k(waven), colden)
        # colden = N_star array.
        # waven = N_wavelengths array
        # tau = N_wave X N_star
        # Flux = N_star X N_wave  ... (may be the other way around.)
        # print("tau", tau)
        # print(np.sum(Flux,axis=1))

        F_ext = [ff * np.power(10, (-1 * tt)) for ff, tt in zip(Flux, tau.T)]


        print(len(F_ext))
        print(len(F_ext[0]))
        Flex_ext = np.zeros(len(F_ext[0]))
        l = 0
        while l < len(F_ext):
            
            Flex_ext += F_ext[l]*star['m'][l]
            l += 1
            print(l/len(F_ext)*100)

        print(Flex_ext)

        Flex = Flex * Lum_sun * 1e3 / np.sum(div)
        Flex_ext = Flex_ext * Lum_sun * 1e3 / np.sum(div)
        waven *= 1e4
        print(waven)


        Vsort = np.where( (waven == 5490) )
        A_lambda = -2.5*np.log10(Flex_ext/Flex)
        A_v = A_lambda[Vsort]
        print(Vsort)
        print(A_lambda)
        print(A_v)
        

        import matplotlib.pyplot as plt
        
        fig1 = plt.figure(figsize=(6,6) )
        f1 = fig1.add_subplot(111)

        #plt.plot(waven, Flex, c='blue')
        #plt.plot(waven, Flex_ext, c='red')
        
        f1.plot(waven,A_lambda/A_v,c='b')
        
        f1.set_xlim(1000,5000)
        f1.set_ylim(0,8)

        f1.set_xlabel("wavelength[Angstrom]")
        f1.set_ylabel("A${_l}{_a}{_m}{_b}{_d}{_a}/A{_V}$")


        fig2 = plt.figure(figsize=(6,6))
        f2 = fig2.add_subplot(111)

        f2.plot(waven, Flex, c='blue')
        f2.plot(waven, Flex_ext, c='red')

        f2.set_xlim(1000,10000)

        f2.set_xlabel("wavelength[Angstrom]")
        f2.set_ylabel("Energy?Flux?[erg]")


        #plt.show()


        # print(np.sum(F_ext, axis = 1))

        # print("before ext", np.sum(Flux, axis=1))
        # print("After ext", np.sum(F_ext, axis=1))
        t6 = time()
        # print("age metal digitize {:.3f}".format(t1-t0))
        # print("setup done {:.3f}".format(t2-t0))
        # print("interpolation coefficient {:.3f}".format(t3-t0))
        # print("flux {:.3f}".format(t4-t0))
        # print("filter, flux, div {:.3f}".format(t5-t0))
        # print("After dust attanuation {:.3f}".format(t6-t0))

        import matplotlib.pyplot as plt
        #print(len(F_ext[0]))
        forsed = []
        i = 0
        while i < len(tau):
            forsed.append(tau[i][0])
            i = i + 1
        # print(waven, tau)
        # if filter_name == 'u':
        # plt.plot(waven, Flux[0], c = 'black')
        # plt.plot(waven, F_ext[0], c = 'red')
        # plt.plot(waven, np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"], c = 'black')
        # plt.scatter(np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"],\

        # plt.legend()
        print("After Extinction = %s" % np.sum((np.sum(F_ext, axis=1) / np.sum(div) * Lum_sun * star['m'])))
        print('')

        return np.sum(Flux, axis=1) / np.sum(div) * Lum_sun * star["m"] , np.sum(F_ext, axis=1) / np.sum(div) * Lum_sun * star["m"]






def ext_curve_k(lam, Rv=4.05):  # calzetti +00 !! #4.05 +- 0.80
    H_frac = 0.76

    lambda1 = 0.48613  # 1e-6 m
    lambda2 = 0.65628  # 1e-6 m

    # Why use only lambda 2
    inv_lambda1 = 1. / lam[lam < lambda2]
    inv_lambda2 = 1. / lam[lam > lambda2]

    k = np.zeros(len(lam))
    # lambda : 0.09 ~ 0.63 1e-6m
    k1 = 2.659 * (-2.156 + (1.509 * inv_lambda1)
                  - (0.198 * (inv_lambda1 ** 2))
                  + (0.011 * (inv_lambda1 ** 3))) + Rv  # Calzetti +94 (25)
    # lambda : 0.63 ~ 5.08 1e-6m
    k2 = 2.659 * (-1.857 + (1.040 * inv_lambda2)) + Rv

    k[lam < lambda2] = k1
    k[lam > lambda2] = k2

    return k




def ext_curve_mw(lam, Rv=3.1):
    Es = 0.44 # EBV

    lambda1 = 1/5.9 # 1e-6 m / may not be used
    lambda2 = 1/3.3 # 1e-6 m 
    lambda3 = 1/1.1 # 1e-6 m
    lambda4 = 1/0.3 # 1e-6 m / may not be used

    lamsort1a = np.where( (lam <= lambda1) )
    lamsort1b = np.where( (lam <= lambda2) & (lam > lambda1) )
    lamsort2 = np.where( (lam > lambda2) & (lam < lambda3) )
    lamsort3 = np.where( (lam > lambda3) )

    inv_lambda1a = 1. / lam[lamsort1a]
    inv_lambda1b = 1. / lam[lamsort1b]
    inv_lambda2 = 1. / lam[lamsort2]
    inv_lambda3 = 1. / lam[lamsort3]

    
    k = np.zeros(len(lam))

    ### UV-a ###
    Fa = -0.04473*( (inv_lambda1a-5.9)**2 ) - 0.009779*( (inv_lambda1a-5.9)**3 )
    Fb = 0.2130*( (inv_lambda1a-5.9)**2 ) +0.1207*( (inv_lambda1a)**3 )
    a1a = 1.752 -0.316*inv_lambda1a -0.104/(((inv_lambda1a-4.67)**2)+0.341) + Fa
    b1a = -3.090 +1.825*inv_lambda1a +1.206*inv_lambda1a/(((inv_lambda1a-4.62)**2)+0.263) + Fb
    k1a = (a1a + b1a/Rv) /Es

    ### UV-b ###
    a1b = 1.752 -0.316*inv_lambda1b -0.104/(((inv_lambda1b-4.67)**2)+0.341) 
    b1b = -3.090 +1.825*inv_lambda1b +1.206*inv_lambda1b/(((inv_lambda1b-4.62)**2)+0.263) 
    k1b = (a1b + b1b/Rv) /Es

    ### Opt & NIR ###
    y = inv_lambda2 -1.82
    a2 = 1 + 0.17699*(y**1) -0.50447*(y**2) -0.02427*(y**3) +0.72085*(y**4) +0.01979*(y**5) -0.77530*(y**6) +0.32999*(y**7)
    b2 = 1.41338*(y**1) +2.28305*(y**2) +1.07233*(y**3) -5.38434*(y**4) -0.62251*(y**5) +5.30260*(y**6) -2.09002*(y**7)
    k2 = (a2 + b2/Rv) /Es

    ### IR ###
    a3 = 0.574*(inv_lambda3 ** 1.61)
    b3 = -0.527*(inv_lambda3 ** 1.61)
    k3 = (a3 + b3/Rv) /Es


    k[lamsort1a] = k1a
    k[lamsort1b] = k1b
    k[lamsort2] = k2
    k[lamsort3] = k3

    return k

#####################################################################



def get_star_colden(nout, idgal, boxlim):
    directory = '/storage1/NewHorizon/'
    gals = load.rd_GM.rd_gal(nout, idgal, wdir=directory)

    s = load.sim.Sim(nout, base=directory)

    gals.header['xg'] = gm2code(gals.header['xg'], gals.info)
    gals.star['x'] = gm2code(gals.star['x'], gals.info)
    gals.star['y'] = gm2code(gals.star['y'], gals.info)
    gals.star['z'] = gm2code(gals.star['z'], gals.info)

    xc, yc, zc = gals.header['xg']

    print('')
    print("Loading star particles...")

    stars = gals.star
    
    xyzsort = np.where( (abs(stars["x"]-xc) < boxlim/gals.info.boxtokpc) &
                (abs(stars["y"]-yc) < boxlim/gals.info.boxtokpc) &
                 (abs(stars["z"]-zc) < boxlim/gals.info.boxtokpc) )
    
    stars = stars[xyzsort]
    print(len(stars))

    stars["x"] -= xc
    stars["y"] -= yc
    stars["z"] -= zc



    stars["vx"] -= gals.header["vg"][0]
    stars["vy"] -= gals.header["vg"][1]
    stars["vz"] -= gals.header["vg"][2]

    theta_xy = 0
    theta_xz = 0
    theta_yz = 0

    # print(star["x"][0])

    stars["x"], stars["y"] = rotation_xy(theta_xy, stars["x"], stars["y"])
    stars["x"], stars["z"] = rotation_xz(theta_xz, stars["x"], stars["z"])
    stars["y"], stars["z"] = rotation_yz(theta_yz, stars["y"], stars["z"])
    stars["vx"], stars["vy"] = rotation_xy(theta_xy, stars["vx"], stars["vy"])
    stars["vx"], stars["vz"] = rotation_xy(theta_xz, stars["vx"], stars["vz"])
    stars["vy"], stars["vz"] = rotation_xy(theta_yz, stars["vy"], stars["vz"])

    ang_x = stars["m"] * (stars["y"] * stars["vz"] - stars["z"] * stars["vy"])
    ang_y = stars["m"] * (stars["z"] * stars["vx"] - stars["x"] * stars["vz"])
    ang_z = stars["m"] * (stars["x"] * stars["vy"] - stars["y"] * stars["vx"])

    params_ang = [np.mean(ang_x),
                  np.mean(ang_y),
                  np.mean(ang_z)]

    stars["x"], stars["y"], stars["z"], theta_xys, theta_yzs = auto_rotation_np(params_ang,
                                                                             stars["x"], stars["y"], stars["z"])

    return theta_xys, theta_yzs




###################################################

def rotation(array1, array2, theta):
    theta *= np.pi / 180
    new_array1 = np.cos(theta) * array1 - np.sin(theta) * array2
    new_array2 = np.sin(theta) * array1 + np.cos(theta) * array2
    return new_array1, new_array2

###############################################

def get_colden(theta_xy, theta_xz, theta_yz, directory=None, file_name=None, quick=False, gridrate=0.5, shift=[0, 0, 0], draw=False, save=False):
    """
    Rotate gas into arbitrary direction
    """

    if gridrate < 2**(-6):
        boxsize=3*10**2
    elif gridrate < 2**(-7):
        boxsize=10**2
    else:
        boxsize = 10**4
    x = np.random.randint(1000, size=boxsize)
    y = np.random.randint(1000, size=boxsize)
    z = np.random.randint(1000, size=boxsize)
    gridsize = 1000 * gridrate  #### notice that gridsize is a half of box's side length
    x, y, z = x - 500, y - 500, z - 500

    x, y = rotation(x, y, theta_xy)
    x, z = rotation(x, z, theta_xz)
    y, z = rotation(y, z, theta_yz)

    x, y, z = x + shift[0], y + shift[1], z + shift[2]

    dsort = np.where((np.sqrt(np.square(x) + np.square(y)) < gridsize * np.sqrt(2))
                     & (abs(x) <= gridsize) & (abs(y) <= gridsize))


    if draw == True:
        plt.show()
    else:
        pass

    z_sort = np.where( (abs(z) <= gridsize) )
    X_zsorted = x[z_sort]
    Y_zsorted = y[z_sort]

    min_xshift = min(X_zsorted)/2/gridsize
    max_xshift = max(X_zsorted)/2/gridsize
    min_yshift = min(Y_zsorted)/2/gridsize
    max_yshift = max(Y_zsorted)/2/gridsize

    min_xshi = 0
    min_yshi = 0
    min_zshi = 0
    max_xshi = 0
    max_yshi = 0
    max_zshi = 0

    if quick == False:

        repeat_size = 20
        i = -repeat_size
        while i <= repeat_size:
            j = -repeat_size
            while j <= repeat_size:
                k = -repeat_size
                checknum = -1
                while k <= repeat_size:

                    if len(np.where((abs(x+2*gridsize*i) <= gridsize) & (abs(y+2*gridsize*j) <= gridsize)& (abs(z+2*gridsize*k)<=gridsize) )[0]) > 0:
                        if checknum == -1:
                            min_xshi = min(min_xshi, i)
                            min_yshi = min(min_yshi, j)
                            min_zshi = min(min_zshi, k)
                            checknum = 0
                            if min_xshi == i or min_yshi ==j or min_zshi ==k:
                                print("min_x_shift is", min_xshi)
                                print("min_y_shift is", min_yshi)
                                print("min_z_shift is", min_zshi,"\n")
                        else:
                            pass
                    else:
                        if checknum == 0:
                            max_xshi = max(max_xshi, i)
                            max_yshi = max(max_yshi, j)
                            max_zshi = max(max_zshi, k)
                            checknum = 1

                            if max_xshi == i or max_yshi == j or max_zshi == k:
                                print("max_x_shift is", max_xshi)
                                print("max_y_shift is", max_yshi)
                                print("max_z_shift is", max_zshi,"\n")
                        else:
                            pass

                    k = k + 1
                j = j +1
            print(100*abs(i+repeat_size)/repeat_size/2)
            i = i +1

            if i == repeat_size+1:
                print("############## Result #################")
                print("min_x_shift is", min_xshi)
                print("min_y_shift is", min_yshi)
                print("min_z_shift is", min_zshi, "\n")

                print("max_x_shift is", max_xshi)
                print("max_y_shift is", max_yshi)
                print("max_z_shift is", max_zshi)
                print("########################################", "\n")

        min_xshi, min_yshi, min_zshi = -1000*np.sqrt(3)/gridsize/2/2,-1000*np.sqrt(3)/gridsize/2/2,-1000*np.sqrt(3)/gridsize/2/2
        max_xshi, max_yshi, max_zshi = 1000*np.sqrt(3)/gridsize/2/2, 1000*np.sqrt(3)/gridsize/2/2, 1000*np.sqrt(3)/gridsize/2/2
    
        print("Only support False/True")
        return

    base_grid_ddx = int(max(max_xshi, abs(min_xshi)))+1
    base_grid_ddy = int(max(max_yshi, abs(min_yshi)))+1
    base_grid_ddz = int(max(max_zshi, abs(min_zshi)))+1
    print("\n","#####################################","\n","base_grid_ddx is ",base_grid_ddx,"\n","#####################################","\n")

    base_grid = np.zeros([2*base_grid_ddz+2+1 ,2*base_grid_ddy+1,2*base_grid_ddx+1])

    i = -base_grid_ddx
    while i <= base_grid_ddx:
        j = -base_grid_ddy
        while j <= base_grid_ddy:
            k = -base_grid_ddz
            while k <= base_grid_ddz:

                component_ijk = len(np.where((abs(x + 2 * gridsize * i) <= gridsize) &
                                    (abs(y + 2 * gridsize * j) <= gridsize) &
                                    (abs(z + 2 * gridsize * k) <= gridsize))[0])/boxsize

                #print(component_ijk)
                base_grid[0][j+base_grid_ddy][i+base_grid_ddx] = i
                base_grid[1][j+base_grid_ddy][i+base_grid_ddx] = j
                base_grid[k+base_grid_ddz+2][j+base_grid_ddy][i+base_grid_ddx] = component_ijk
                #base_grid[i+base_grid_ddx][j+base_grid_ddy][k+base_grid_ddz] = component_ijk
                k = k + 1
            j = j +1
        print(100*abs(i+base_grid_ddx)/base_grid_ddx/2)
        i = i +1

    print(base_grid)

    if save == True:
        save_route = directory#"/home/jangjk/PycharmProjects/PP/week3/base_grid/"
        #file_name = "first_try_basegrid.npy"
        route_name = save_route+file_name

        np.save(route_name,base_grid)

    return len(dsort[0]), base_grid  # /grid/grid/1000


def smoothing(theta_xy,theta_xz,theta_yz, directory, file_name, quick, gridrate=0.5, shift=[0, 0, 0], plt=False, save=False):
    A = []
    Base_grid = []
    i = 0
    while i < len(theta_xz):
        if i == 0:
            a, ith_base_grid = get_colden(theta_xy, theta_xz[i], theta_yz, gridrate=gridrate, shift=shift, draw=plt, quick=quick,
                                          directory =directory,file_name=file_name, save=save )
        else:
            a, ith_base_grid = get_colden(theta_xy, theta_xz[i], theta_yz, gridrate=gridrate, shift=shift, quick=quick,
                                          directory=directory, file_name=file_name ,save=save)
        A.append(a)

        if len(theta_xz) == 1:
            Base_grid = ith_base_grid
        elif len(theta_xz) > 1:
            Base_grid.append(ith_base_grid)
        else:
            pass

        i += 1
        print(np.round(100 * i / len(theta_xz), 2))

    return np.array(A), Base_grid



###############################################################

def load_gal(nout, idgal,theta_xz_1, theta_yz_1, save_directory,fixed_idgal,
             boxlim, celllim, boxlim_xbot, boxlim_xupp, boxlim_ybot, boxlim_yupp,
             drawlim_xbot, drawlim_xupp, drawlim_ybot, drawlim_yupp,drawpix,
             saving=True, fuv=False, BVR=False, SDSS=False, JHK=False):

    directory = '/storage1/NewHorizon/'
    directory5 = '/storage5/NewHorizon/'
    directory_galactica = '/storage5/Galactica/galaxy_14667/TREE_STARS/'
    try:
        gal = load.rd_GM.rd_gal(nout, idgal, wdir=directory)
        s = load.sim.Sim(nout, base=directory)
    except:
        try:    
            gal = load.rd_GM.rd_gal(nout, idgal, wdir=directory5)
            s = load.sim.Sim(nout, base=directory5)
        except:
            gal = load.rd_GM.rd_gal(nout, idgal, wdir=directory_galactica,fname="galactica")
            s = load.sim.Sim(nout, base=directory_galactica,data_dir='galactica')
            

    gal.header['xg'] = gm2code(gal.header['xg'], gal.info)
    gal.star['x'] = gm2code(gal.star['x'], gal.info)
    gal.star['y'] = gm2code(gal.star['y'], gal.info)
    gal.star['z'] = gm2code(gal.star['z'], gal.info)

    xc, yc, zc = gal.header['xg']

    star = gal.star


    xyzsort = np.where( 
                        (abs(star["x"]-xc) <  boxlim/gal.info.boxtokpc) & 
                        (abs(star["y"]-yc) <  boxlim/gal.info.boxtokpc) & 
                        (abs(star["z"]-zc) <  boxlim/gal.info.boxtokpc)  
                      )

    print('')
    print("Loading star particles...")


    star = star[xyzsort]

    star["x"] -= xc
    star["y"] -= yc
    star["z"] -= zc


    star['vx'] -= gal.header['vg'][0]
    star['vy'] -= gal.header['vg'][1]
    star['vz'] -= gal.header['vg'][2]


    ###############################################################################################3

    # radius = 0.5 * max([gal.star['x'].ptp(), gal.star['y'].ptp(), gal.star['z'].ptp()])
    radius = celllim / gal.info.boxtokpc
    Rlim_sim = radius  # region setting


    s.set_ranges([[xc - Rlim_sim, xc + Rlim_sim], [yc - Rlim_sim, yc + Rlim_sim], [zc - Rlim_sim, zc + Rlim_sim]])


    ############################################################################################
    # Load stellar particle
    s.add_part(ptypes=['dm id pos vel mass'],fortran=True)

    dm = s.part.dm
    dmsort = np.where( (dm['m']>0) & (dm['id'] >=0) )
    dm = dm[dmsort]

    
    dm["x"] -= xc
    dm["y"] -= yc
    dm["z"] -= zc
    dm["m"] *= s.info.msun

    #dm["x"] *= s.info.boxtokpc
    #dm["y"] *= s.info.boxtokpc
    #dm["z"] *= s.info.boxtokpc


    ################################################################################################
    # Load cell data
    s.add_hydro()

    #################### load gas cells #############################
    print('')
    print('Loading gas cells...')

    ind_g = np.where(np.square(s.hydro.cell['x'] - xc) +
                     np.square(s.hydro.cell['y'] - yc) +
                     np.square(s.hydro.cell['z'] - zc) < np.square(Rlim_sim))[0]

    cell = s.hydro.cell[ind_g]
    #print(cell["x"])
    #print(len(cell["x"]))
    
    cell['x'] -= xc
    cell['y'] -= yc
    cell['z'] -= zc


    cell['var1'] = cell['var1'] * s.info.kms - gal.header['vg'][0]
    cell['var2'] = cell['var2'] * s.info.kms - gal.header['vg'][1]
    cell['var3'] = cell['var3'] * s.info.kms - gal.header['vg'][2]



    cell_rho = cell['var0'] * s.info.unit_nH
    cell_T = cell['var4'] / cell['var0'] * s.info.unit_T2
    cell_mgas = cell['var0'] * s.info.unit_d * (cell['dx'] * s.info.boxtokpc * kpc_in_cm) ** 3 / msun_in_g

    print('complete')
    

    #######################################################################
    # Calcuate Nvec direction
    print('')
    print("Loading star particles...")


    ang_x = star["m"] * (star["y"] * star["vz"] - star["z"] * star["vy"])
    ang_y = star["m"] * (star["z"] * star["vx"] - star["x"] * star["vz"])
    ang_z = star["m"] * (star["x"] * star["vy"] - star["y"] * star["vx"])

    params_ang = [np.mean(ang_x),
                  np.mean(ang_y),
                  np.mean(ang_z)]



    star["x"], star["y"], star["z"], theta_xy_0,theta_yz_0 = auto_rotation_np(params_ang,
                                                           star["x"], star["y"], star["z"],thetype=1)

    cell["x"], cell["y"], cell["z"], theta_xy_0, theta_yz_0 = auto_rotation_np(params_ang,
                                                           cell["x"], cell["y"], cell["z"],thetype=1)

    cell["var1"], cell["var2"], cell["var3"], theta_xy_0, theta_yz_0 = auto_rotation_np(params_ang, cell["var1"], cell["var2"],
                                                                    cell["var3"],thetype=1)

    CTC = np.vectorize(cartesian_to_cylinder)

    print('')
    print("Initiate Rotation")

    x, y, z = star['x'],star['y'],star['z']
    vx, vy, vz = auto_rotation_np(params_ang, star["vx"], star["vy"], star["vz"])
    r, phi, z, v_r, v_phi, v_z = CTC(x, y, z, vx, vy, vz)

    print("Star particles complete")
    print("")

    dm_x, dm_y, dm_z = auto_rotation_np(params_ang, dm["x"], dm["y"], dm["z"])
    dm_vx, dm_vy, dm_vz = auto_rotation_np(params_ang, dm["vx"], dm["vy"], dm["vz"])
    dm_r, dm_phi, dm_z, dm_v_r, dm_v_phi, dm_v_z = CTC(dm_x, dm_y, dm_z, dm_vx, dm_vy, dm_vz)

    print("DM particles complete")
    print("")

    cell_x, cell_y, cell_z = cell["x"], cell["y"], cell["z"]
    cell_vx, cell_vy, cell_vz = cell["var1"], cell["var2"], cell["var3"]
    cell_r, cell_phi, cell_z, cell_v_r, cell_v_phi, cell_v_z = CTC(cell_x, cell_y, cell_z, cell_vx, cell_vy, cell_vz)

    print("Gas Cells complete")
    print("")


    ##################################################

    star_mass = star['m']*1e11*msun_in_g
    dm_mass = dm['m']*msun_in_g
    gas_mass = cell_mgas * msun_in_g

    x *= gal.info.boxtokpc
    y *= gal.info.boxtokpc
    z *= gal.info.boxtokpc
    r *= gal.info.boxtokpc
    
    cell_x *= gal.info.boxtokpc
    cell_y *= gal.info.boxtokpc
    cell_z *= gal.info.boxtokpc
    cell_r *= gal.info.boxtokpc

    dm_x *= gal.info.boxtokpc
    dm_y *= gal.info.boxtokpc
    dm_z *= gal.info.boxtokpc
    dm_r *= gal.info.boxtokpc

    x *= kpc_in_cm
    y *= kpc_in_cm
    z *= kpc_in_cm
    r *= kpc_in_cm

    dm_x *= kpc_in_cm
    dm_y *= kpc_in_cm
    dm_z *= kpc_in_cm
    dm_r *= kpc_in_cm

    cell_x *= kpc_in_cm
    cell_y *= kpc_in_cm
    cell_z *= kpc_in_cm
    cell_r *= kpc_in_cm
    
    vx *= 1e5
    vy *= 1e5
    vz *= 1e5

    v_r *= 1e5
    v_phi *= 1e5
    v_z *= 1e5

    ##################################################
    T_seperate = 6+0.25*np.log10(cell_rho)

    cold_gas_sort = np.where(  (np.log10(cell_T) < T_seperate)   )

    cold_gas_mass = cell_mgas[cold_gas_sort]


    ang_cell = cell_r * cell_v_phi * cell_mgas
    ang_cell = np.log10(ang_cell)


    plt.hist(ang_cell[cold_gas_sort],bins=500,histtype='step',color='b')
    plt.hist(ang_cell,bins=500,histtype='step',color='orange')
    #plt.show()







    sig_r = ( sum( ( v_r - sum(v_r)/len(v_r) )**2 ) / len(v_r) )**0.5
    sig_phi = ( sum( ( v_phi - sum(v_phi)/len(v_phi) )**2 ) / len(v_phi) )**0.5
    sig_z = ( sum( ( v_z - sum(v_z)/len(v_z) )**2 ) / len(v_z) )**0.5

    sig_r = np.std(v_r)
    sig_phi = np.std(v_phi)
    sig_z = np.std(v_z)
    sig_cylinder =  ( (sig_r**2 + sig_phi**2 + sig_z**2) /3 )**0.5

    VoS = abs(np.mean(v_phi)) / sig_cylinder

    print("")
    print( "V/sigma = %s" %(VoS) )
    print('')

    ##################################################
    tc = utils.cosmology.Timeconvert(info=gal.info)
    starage = tc.time2gyr(star['time'],z_now = gal.info.zred)

    r_new = list(r) + list(dm_r) + list(cell_r)
    z_new = list(z) + list(dm_z) + list(cell_z)
    mass_new = list(star_mass) + list(dm_mass) + list(gas_mass)

    r_new = np.array(r_new)
    z_new = np.array(z_new)
    mass_new = np.array(mass_new)

    radsort = np.argsort(np.square(r_new)+np.square(z_new))
    r_new = r_new[radsort]
    z_new = z_new[radsort]
    mass_new2 = mass_new[radsort]

    cum_mass = np.cumsum(mass_new2)

    dist = np.sqrt(np.square(r)+np.square(z))
    dist_bin=np.sqrt(np.square(r_new)+np.square(z_new))

    rad_index = np.digitize(dist,bins=dist_bin)-1
    rsort = np.where( (r>0) )

    KE = 0.5*star_mass[rsort]*(np.square(v_phi[rsort])+np.square(v_z[rsort])+np.square(v_r[rsort]))
    PE = -1*G*cum_mass[rad_index][rsort]*star_mass[rsort]/np.sqrt(np.square(r[rsort])+np.square(z[rsort]))
    Etot = KE + PE

    v_cir = np.sqrt(-2*Etot/star_mass[rsort])

    J_z = r[rsort]*v_phi[rsort]
    J_cir = G*cum_mass[rad_index][rsort]/v_cir

    cir_param = J_z/J_cir

    print(cir_param)

    disc_sort = np.where( (cir_param > 0.5) )
    sph_sort = np.where( (cir_param < 0.5) )


    DtT = len(disc_sort[0])/(len(sph_sort[0])+len(disc_sort[0]))
    print(DtT)

    plt.scatter(r[rsort]/kpc_in_cm,cir_param,alpha=0.5,c=starage[rsort],cmap='jet')
    #plt.savefig("./image/%s/%s" %(fixed_idgal,str(nout)+".png"))

    ##################################################

    x /= gal.info.boxtokpc
    y /= gal.info.boxtokpc
    z /= gal.info.boxtokpc
    r /= gal.info.boxtokpc

    cell_x /= gal.info.boxtokpc
    cell_y /= gal.info.boxtokpc
    cell_z /= gal.info.boxtokpc
    cell_r /= gal.info.boxtokpc

    dm_x /= gal.info.boxtokpc
    dm_y /= gal.info.boxtokpc
    dm_z /= gal.info.boxtokpc
    dm_r /= gal.info.boxtokpc

    x /= kpc_in_cm
    y /= kpc_in_cm
    z /= kpc_in_cm

    dm_x /= kpc_in_cm
    dm_y /= kpc_in_cm
    dm_z /= kpc_in_cm

    cell_x /= kpc_in_cm
    cell_y /= kpc_in_cm
    cell_z /= kpc_in_cm

    vx /= 1e5
    vy /= 1e5
    vz /= 1e5





    ##################################################

    xyzsort = np.where( ( star["x"] >= boxlim_xbot/gal.info.boxtokpc ) &
                        ( star["x"] <  boxlim_xupp/gal.info.boxtokpc ) &
                        ( star["y"] >= boxlim_ybot/gal.info.boxtokpc ) &
                        ( star["y"] <  boxlim_yupp/gal.info.boxtokpc ) )

            
    star = star[xyzsort]




    ###################################################

    fixed_grid_size = cell_ddx_grid_size[0] # dx_min

    changed_grid_size0 = cell_ddx_grid_size[0]
    changed_grid_size1 = cell_ddx_grid_size[1]
    changed_grid_size2 = cell_ddx_grid_size[2]
    changed_grid_size3 = cell_ddx_grid_size[3]
    changed_grid_size4 = cell_ddx_grid_size[4]
    changed_grid_size5 = cell_ddx_grid_size[5]
    changed_grid_size6 = cell_ddx_grid_size[6]
    changed_grid_size7 = cell_ddx_grid_size[7]

    grid_rate0 = fixed_grid_size / changed_grid_size0 / 2
    grid_rate1 = fixed_grid_size / changed_grid_size1 / 2
    grid_rate2 = fixed_grid_size / changed_grid_size2 / 2
    grid_rate3 = fixed_grid_size / changed_grid_size3 / 2
    grid_rate4 = fixed_grid_size / changed_grid_size4 / 2
    grid_rate5 = fixed_grid_size / changed_grid_size5 / 2
    grid_rate6 = fixed_grid_size / changed_grid_size6 / 2
    grid_rate7 = fixed_grid_size / changed_grid_size7 / 2

    grid_size0 = int(1000 * grid_rate0)
    grid_size1 = int(1000 * grid_rate1)
    grid_size2 = int(1000 * grid_rate2)
    grid_size3 = int(1000 * grid_rate3)
    grid_size4 = int(1000 * grid_rate4)
    grid_size5 = int(1000 * grid_rate5)
    grid_size6 = int(1000 * grid_rate6)
    grid_size7 = int(1000 * grid_rate7)

    shift0 = [0, 0, 0]
    shift1 = [0, 0, 0]
    shift2 = [0, 0, 0]
    shift3 = [0, 0, 0]
    shift4 = [0, 0, 0]
    shift5 = [0, 0, 0]
    shift6 = [0, 0, 0]
    shift7 = [0, 0, 0]

    #####################################################

    size = 0
    plts = True
    plts = False
    quick = True
    save = False
    directories = "/home/jangjk816/Project/Mock/week7/base_grid/"

    ######################################################

    O, rot_cell0 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate0, shift=shift0, plt=plts, quick=quick,
                             directory=directories, file_name="0_basegrid_quick.npy", save=save)

    A, rot_cell1 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate1, shift=shift1, plt=plts, quick=quick,
                             directory=directories, file_name="1_basegrid_quick.npy", save=save)

    B, rot_cell2 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate2, shift=shift2, plt=plts, quick=quick,
                             directory=directories, file_name="2_basegrid_quick.npy", save=save)

    C, rot_cell3 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate3, shift=shift3, plt=plts, quick=quick,
                             directory=directories, file_name="3_basegrid_quick.npy", save=save)

    D, rot_cell4 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate4, shift=shift4, plt=plts, quick=quick,
                             directory=directories, file_name="4_basegrid_quick.npy", save=save)

    E, rot_cell5 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate5, shift=shift5, plt=plts, quick=quick,
                             directory=directories, file_name="5_basegrid_quick.npy", save=save)


    F, rot_cell6 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                                         gridrate=grid_rate6, shift=shift6, plt=plts, quick=quick,
                                                                      directory=directories, file_name="6_basegrid_quick.npy", save=save)

    Gs, rot_cell7 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                                         gridrate=grid_rate7, shift=shift7, plt=plts, quick=quick,
                                                                      directory=directories, file_name="7_basegrid_quick.npy", save=save)



    ######################################################################
    theta_xz = theta_xz_1[0]
    star["x"], star["z"] = rotation(star["x"], star["z"],theta_xz)
    cell["x"], cell["z"] = rotation(cell["x"], cell["z"],theta_xz)

    
    theta_yz = theta_yz_1[0]
    star["y"], star["z"] = rotation(star["y"], star["z"],theta_yz)
    cell["y"], cell["z"] = rotation(cell["y"], cell["z"],theta_yz)



    print("\n","rotation complete","\n")

    ########################################################################


    cell_dx_min = np.min(cell["dx"])
    cell_dx_max = np.max(cell["dx"])
    
    ddx_h = cell["dx"] * 0.5
    print("min_ddx_h =",min(ddx_h),"\n","max_ddx_h =",max(ddx_h))

    dx1sort = np.where((ddx_h >= 1e-8) & (ddx_h < 1e-7))
    dx2sort = np.where((ddx_h >= 1e-7) & (ddx_h < 1e-6))
    dx3sort = np.where((ddx_h >= 1e-6) & (ddx_h < 1e-5))

    print("Number of cells in ddx:[1e-8,1e-7) = %s" % (len(ddx_h[dx1sort])))
    print("Number of cells in ddx:[1e-7,1e-6) = %s" % (len(ddx_h[dx2sort])))
    print("Number of cells in ddx:[1e-6,1e-5) = %s" % (len(ddx_h[dx3sort])))
    print("total number of cell = %s" % (len(cell)))
    print("total number of selected cell = %s" % (len(ddx_h[dx1sort]) + len(ddx_h[dx2sort]) + len(ddx_h[dx3sort])))    


    sub_cell = cell[((cell["x"] + ddx_h) > star["x"].min()) * (cell["x"] - ddx_h < star["x"].max()) * \
                    ((cell["y"] + ddx_h) > star["y"].min()) * (cell["y"] - ddx_h < star["y"].max()) * \
                    ((cell["z"] + ddx_h) < star["z"].max())]

    sub_cell_T = sub_cell['var4'] / sub_cell['var0'] * s.info.unit_T2


    print("len subcell", len(sub_cell))
    if len(sub_cell) < 1:
        print("something's wrong with sub_cell")
        return 


    ########################################################################


    sub_ddx_h = sub_cell["dx"] * 0.5
    sub_ddx_list = np.unique(sub_ddx_h)
    print("\n","sub_cell's dx/2 size =",sub_ddx_list,"\n")



    xl = sub_cell["x"] - sub_ddx_h*np.sqrt(3)
    xr = sub_cell["x"] + sub_ddx_h*np.sqrt(3)
    yl = sub_cell["y"] - sub_ddx_h*np.sqrt(3)
    yr = sub_cell["y"] + sub_ddx_h*np.sqrt(3)
    zl = sub_cell["z"] - sub_ddx_h*np.sqrt(3)
    zr = sub_cell["z"] + sub_ddx_h*np.sqrt(3)

    print("min_xl = %s, max_xl = %s" % (min(xl), max(xl)))
    print("min_xr = %s, max_xr = %s" % (min(xr), max(xr)))
    print("min_yl = %s, max_yl = %s" % (min(yl), max(yl)))
    print("min_yr = %s, max_yr = %s" % (min(yr), max(yr)))
    print("min_zl = %s, max_zl = %s" % (min(zl), max(zl)))
    print("min_zr = %s, max_zr = %s" % (min(zr), max(zr)))

    xrange = (xl.min(), xr.max())
    yrange = (yl.min(), yr.max())
    zrange = (zl.min(), zr.max())
    xspan = xrange[1] - xrange[0]
    yspan = yrange[1] - yrange[0]
    zspan = zrange[1] - zrange[0]

    ############################################################
    vmin_x, vmax_x = xrange[0], xrange[1]
    vmin_y, vmax_y = yrange[0], yrange[1]
    vmin_z, vmax_z = zrange[0], zrange[1]
    ##############################################################

    ssort = np.where((star["x"] >= vmin_x) & (star["x"] <= vmax_x) &
                     (star["y"] >= vmin_y) & (star["y"] <= vmax_y))

    
    print("star_x_min and max boxsize =", min(star["x"][ssort]), max(star["x"][ssort]))
    print("star_y_min and max boxsize =", min(star["y"][ssort]), max(star["y"][ssort]))

    x_real_span = max(star["x"][ssort]) - min(star["x"][ssort])
    y_real_span = max(star["y"][ssort]) - min(star["y"][ssort])
    z_real_span = max(star["z"][ssort]) - min(star["z"][ssort])
    
    x_real_span = x_real_span*gal.info.boxtokpc
    y_real_span = y_real_span*gal.info.boxtokpc
    z_real_span = z_real_span*gal.info.boxtokpc

    print("\n","x_real_span is ",x_real_span)
    print("\n","y_real_span is ",y_real_span)
    print("\n","z_real_span is ",z_real_span)

    print(xrange)
    print(yrange)

    print("\n","xspan is ",xspan*gal.info.boxtokpc)
    print("\n","yspan is ",yspan*gal.info.boxtokpc) 
    print("\n","zspan is ",zspan*gal.info.boxtokpc)

    npixx = np.int(np.ceil(xspan / cell_dx_min))
    npixy = np.int(np.ceil(yspan / cell_dx_min))
    npixz = np.int(np.ceil(zspan / cell_dx_min))
    print("\n","(npixx,npixy,npixz) =","(",npixx,",",npixy,",",npixz,")","\n")


    h = np.histogram2d(star["x"], star["y"],
                       bins=[npixx, npixy],
                       range=[xrange, yrange])

    dxmap = dymap = cell_dx_min

    # Sort stars
    sx = np.searchsorted(h[1], star["x"]) - 1
    sy = np.searchsorted(h[2], star["y"]) - 1

    cx = np.searchsorted(h[1], sub_cell["x"]) - 1
    cy = np.searchsorted(h[2], sub_cell["y"]) - 1

    zindex = np.linspace(zrange[0],zrange[1],npixz)
#    print(zindex)
    
    sz = np.searchsorted(zindex, star["z"]) -1
    cz = np.searchsorted(zindex, sub_cell["z"]) -1
#    print(star["z"])
#    print(zindex[sz])


    rot_cell_length = [len(rot_cell0)-2,len(rot_cell1)-2,
                        len(rot_cell2)-2,len(rot_cell3)-2,
                          len(rot_cell4)-2,len(rot_cell5)-2,
                            len(rot_cell6)-2,len(rot_cell7)-2]
    
    
    
    print(rot_cell_length)
    
    x_padding_layer = int(max(rot_cell_length)/2)
    y_padding_layer = int(max(rot_cell_length)/2)
    z_padding_layer = int(max(rot_cell_length)/2) 
    
    print("... Preparing for stacking column density of cells ...")
    
    Grid = np.zeros([npixz+z_padding_layer*2+2,npixy+y_padding_layer*2,npixx+x_padding_layer*2]) # New gas distribution 
    
    ######### Notice that the order of index is not [x,y,z] but [z,y,x] ###################
    ######### Additional 2 arrays in z component is for x,y indexing when we draw maps ############ 
    
    print("x_length of empty grid = ",len(Grid[0][0]))
    print("y_length of empty grid = ",len(Grid[0]))
    print("z_length of empty grid = ",len(Grid))

    rot_cell_length = np.array(rot_cell_length)
    rot_cell_centre_num = (rot_cell_length-1)/2-1+1+2 ### -1 : indexing is start from zero, not one.
                                                        ### +2 : first and second layer is xy indexing layer.
                                                          ### +1 : length + center(1) + length
    sc = 0 #+len(sub_cell)
    while sc < len(sub_cell):
        rot_cell_index_sc = np.searchsorted(sub_ddx_list,sub_ddx_h[sc])
        if rot_cell_index_sc >= 5:
            pass

        rot_cell_length_sc = rot_cell_length[rot_cell_index_sc]
        rot_cell_centre_num_sc = int(rot_cell_centre_num[rot_cell_index_sc])
        #print(sub_ddx_h[sc],sub_ddx_list[np.searchsorted(sub_ddx_list,sub_ddx_h[sc])])
        rot_cell_sc = []
        if rot_cell_index_sc == 0:
            rot_cell_sc = rot_cell0
        elif rot_cell_index_sc == 1:
            rot_cell_sc = rot_cell1
        elif rot_cell_index_sc == 2:
            rot_cell_sc = rot_cell2
        elif rot_cell_index_sc == 3:
            rot_cell_sc = rot_cell3
        elif rot_cell_index_sc == 4:
            rot_cell_sc = rot_cell4
        elif rot_cell_index_sc == 5:
            rot_cell_sc = rot_cell5
        elif rot_cell_index_sc == 6:
            rot_cell_sc = rot_cell6
        elif rot_cell_index_sc == 7:
            rot_cell_sc = rot_cell7
        else:
            print("Error occured at calculating part of indexing sc'th sub_cell's index")
            break
        
        Gx_sc = cx[sc]+x_padding_layer
        Gy_sc = cy[sc]+y_padding_layer
        Gz_sc = cz[sc]+z_padding_layer+2 ### +2 for removing xy-indexing layer
        
        test_rcsc = np.zeros([2+rot_cell_length_sc,rot_cell_length_sc,rot_cell_length_sc])
        

        rcsc_low = rot_cell_centre_num_sc-int((rot_cell_length_sc-1)/2)
        rcsc_upp = rot_cell_centre_num_sc+int((rot_cell_length_sc-1)/2)

        test_rcsc[rcsc_low:rcsc_upp+1] = rot_cell_sc[rcsc_low:rcsc_upp+1]

        Gx_sc_low = Gx_sc - int((rot_cell_length_sc-1)/2)
        Gx_sc_upp = Gx_sc + int((rot_cell_length_sc-1)/2)

        Gy_sc_low = Gy_sc - int((rot_cell_length_sc-1)/2)
        Gy_sc_upp = Gy_sc + int((rot_cell_length_sc-1)/2)
        
        Gz_sc_low = Gz_sc - int((rot_cell_length_sc-1)/2)
        Gz_sc_upp = Gz_sc + int((rot_cell_length_sc-1)/2)
       

        
        if sub_cell_T[sc] <= 8000:
            Grid[Gz_sc_low:Gz_sc_upp+1,Gy_sc_low:Gy_sc_upp+1,Gx_sc_low:Gx_sc_upp+1] += rot_cell_sc[rcsc_low:rcsc_upp+1]*sub_cell["var0"][sc]*(sub_cell["dx"][sc]**3)/(cell_dx_min**2)

        else:
            pass

        sc += 1
        print( np.round(100*sc/len(sub_cell),2) )
    
    j = 0  
    while j < len(Grid[0]):
        i = 0
        while i < len(Grid[0][0]):
            Grid[0][j][i] = i
            Grid[1][j][i] = j
            i += 1
        j += 1
    
    
    sumup = 0
    k = 2
    while k < len(Grid):
        sumup_k = sum(sum(Grid[k]))
        sumup = sumup + sumup_k
        k = k +1

    file_name_grid = str(nout) + '_' +str(idgal) + '_' + str(theta_xz_1[0])+ '_' + str(theta_yz_1[0]) + '_test_grid.npy'

    #np.save("./base_grid/%s" %(file_name_grid),Grid)
    #print("Grid Saved")
    
    #Grid = np.load("./base_grid/%s" %(file_name_grid))
    ############################## Calculate colden for every star ###########################
    nstar = len(star)
    colden_star = np.zeros(nstar)
    
    i = 0 #+ len(star)
    while i < len(star):
        sx_i = sx[i]
        sy_i = sy[i]
        sz_i = sz[i]
        
        #### 이 부분이 중요함. 방향이 어디이냐에 따라서 나오는 colden 결과가 완전 반대로 나올 수도 있다. ###
        colden_star[i] = sum(Grid[z_padding_layer+2:z_padding_layer+2+sz_i,y_padding_layer+sy_i,x_padding_layer+sx_i])
        #colden_star[i] = sum(Grid[z_padding_layer+2+sz_i+1:len(Grid),y_padding_layer+sy_i,x_padding_layer+sx_i])
        
        print( np.round(100*i/len(star),2) )
        i += 1
    
    colden = colden_star * gal.info.unit_nH * gal.info.unit_l #/ gal.info.boxtokpc
    
    print('')
    print('snapshot Number = %s' %(nout) )
    print('galaxy Number = %s' %(idgal) )
    print('Rotate angle along y-axis = %s' %(theta_xz_1[0]) )
    print('Rotate angle along x-axis = %s' %(theta_yz_1[0]) )
    print("Sum up all Grid componets is ",sumup)
    print('')

    file_name_colden = str(nout) + '_' +str(idgal) + '_' + str(theta_xz_1[0])+ '_' + str(theta_yz_1[0]) + '_temporary_colden.npy'
    #np.save("./colden/%s" %(file_name_colden),colden)
    #print("Colden Saved")
    
    #colden = np.load("./colden/%s" %(file_name_colden))
    # 여기까지!
    #########################################################################
    # Generating Mock image

    #sed_ori = qmc_ori.Simplemock()
    #u_flux = sed_ori.get_flux(gal.star,filter_name='u', info=gal.info)
    #g_flux = sed_ori.get_flux(gal.star,filter_name='g', info=gal.info)
    #r_flux = sed_ori.get_flux(gal.star, filter_name='r', info=gal.info)
    #i_flux = sed_ori.get_flux(gal.star, filter_name='i', info=gal.info)

    maxbin = boxlim/gal.info.boxtokpc

    drawlim_xbot = drawlim_xbot/gal.info.boxtokpc
    drawlim_xupp = drawlim_xupp/gal.info.boxtokpc
    drawlim_ybot = drawlim_ybot/gal.info.boxtokpc
    drawlim_yupp = drawlim_yupp/gal.info.boxtokpc

    npix = drawpix#1000
    weight_R = 1.0
    weight_G = 1.0
    weight_B = 1.0


    ssort = np.where((star["x"] == star["x"]) )
    sed = Simplemock()

    #FUV_flux, FUV_flux_extinc = sed.get_SED(star, num=1, cell=cell, colden=colden, extinction=True,
     #                                        info=gal.info)

    ## FUV ##
    if saving == True and fuv == True:
        sed.__init__(filter_system="FUV",info=gal.info)
        FUV_flux, FUV_flux_extinc = sed.get_flux(star, num=1, filter_name="FUV", cell=cell, colden=colden, extinction = True, info=gal.info)

        FUV = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(),
                         weights=FUV_flux_extinc[ssort], cmap='gray',
                         range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])

        FUV_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(),
                         weights=FUV_flux[ssort], cmap='gray',
                         range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])


        hduFUV = fits.PrimaryHDU(FUV[0] * weight_B)
        hdulFUV = fits.HDUList([hduFUV])
        file_name_FUV = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_FUVfilter.fits'
        hdulFUV.writeto('%s/%s/%s' % (str(save_directory),str(nout), file_name_FUV))

        hduFUV_non = fits.PrimaryHDU(FUV_non[0] * weight_B)
        hdulFUV_non = fits.HDUList([hduFUV_non])
        file_name_FUV_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_FUVfilter_non.fits'
        hdulFUV_non.writeto('%s/%s/%s' % (str(save_directory),str(nout), file_name_FUV_non))

    else:
        pass


    ### Jhonson & Cousin ###
    if saving == True and BVR == True:
        print("what")
        sed.__init__(filter_system="BVR",info=gal.info)

        B_flux, B_flux_extinc = sed.get_flux(star, num=1, filter_name="B", cell=cell, colden=colden, extinction = True, info=gal.info)
        B = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=B_flux_extinc[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        B_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=B_flux[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        
        hduB = fits.PrimaryHDU(B[0] * weight_B  )
        hdulB = fits.HDUList([hduB])
        file_name_B = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_Bfilter.fits'
        hdulB.writeto('%s/%s' % (str(save_directory),file_name_B))

        hduB_non = fits.PrimaryHDU(B_non[0] * weight_B  )
        hdulB_non = fits.HDUList([hduB_non])
        file_name_B_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_Bfilter_non.fits'
        hdulB_non.writeto('%s/%s' % (str(save_directory),file_name_B_non))



        V_flux, V_flux_extinc = sed.get_flux(star, num=1, filter_name="V", cell=cell, colden=colden, extinction = True, info=gal.info)
        V = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=V_flux_extinc[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        V_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=V_flux[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        
        hduV = fits.PrimaryHDU(V[0] * weight_G  )
        hdulV = fits.HDUList([hduV])
        file_name_V = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_Vfilter.fits'
        hdulV.writeto('%s/%s' % (str(save_directory),file_name_V))

        hduV_non = fits.PrimaryHDU(V_non[0] * weight_G  )
        hdulV_non = fits.HDUList([hduV_non])
        file_name_V_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_Vfilter_non.fits'
        hdulV_non.writeto('%s/%s' % (str(save_directory),file_name_V_non))

        
        
        R_flux, R_flux_extinc = sed.get_flux(star, num=1, filter_name="R", cell=cell, colden=colden, extinction = True, info=gal.info)
        R = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=R_flux_extinc[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        R_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=R_flux[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])

        hduR = fits.PrimaryHDU(R[0] * weight_R  )
        hdulR = fits.HDUList([hduR])
        file_name_R = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_Rfilter.fits'
        hdulR.writeto('%s/%s' % (str(save_directory),file_name_R))

        hduR_non = fits.PrimaryHDU(R_non[0] * weight_R  )
        hdulR_non = fits.HDUList([hduR_non])
        file_name_R_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_Rfilter_non.fits'
        hdulR_non.writeto('%s/%s' % (str(save_directory),file_name_R_non))
    
    else:
        pass



    if saving == True and SDSS == True:
        ### SDSS ####
        sed.__init__(filter_system="SDSS",info=gal.info)

        u_flux, u_flux_extinc = sed.get_flux(star, num=1, filter_name='u', cell=cell, colden=colden, extinction = True, info=gal.info)

        u = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(),
                       weights=u_flux_extinc[ssort], cmap='gray',
                       range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        u_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(),
                       weights=u_flux[ssort], cmap='gray',
                       range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])

        hduu = fits.PrimaryHDU(u[0] * weight_B)
        hdulu = fits.HDUList([hduu])
        file_name_u = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(
            theta_yz_1[0]) + '_ufilter.fits'
        hdulu.writeto('%s/%s/%s' % (str(save_directory),str(nout), file_name_u))

        hduu_non = fits.PrimaryHDU(u_non[0] * weight_B)
        hdulu_non = fits.HDUList([hduu_non])
        file_name_u_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(
            theta_yz_1[0]) + '_ufilter_non.fits'
        hdulu_non.writeto('%s/%s/%s' % (str(save_directory),str(nout),file_name_u_non))



        g_flux, g_flux_extinc = sed.get_flux(star, num=1, filter_name='g', cell=cell, colden=colden, extinction = True, info=gal.info)

        g = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(),
                       weights=g_flux_extinc[ssort], cmap='gray',
                       range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        g_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(),
                       weights=g_flux[ssort], cmap='gray',
                       range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])

        hdug = fits.PrimaryHDU(g[0] * weight_B)
        hdulg = fits.HDUList([hdug])
        file_name_g = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(
            theta_yz_1[0]) + '_gfilter.fits'
        hdulg.writeto('%s/%s/%s' % (str(save_directory),str(nout), file_name_g))

        hdug_non = fits.PrimaryHDU(g_non[0] * weight_B)
        hdulg_non = fits.HDUList([hdug_non])
        file_name_g_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(
            theta_yz_1[0]) + '_gfilter_non.fits'
        hdulg_non.writeto('%s/%s/%s' % (str(save_directory),str(nout), file_name_g_non))


        r_flux, r_flux_extinc = sed.get_flux(star, num=1, filter_name='r', cell=cell,colden=colden, extinction = True, info=gal.info)

        r = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(),
                       weights=r_flux_extinc[ssort], cmap='gray',
                       range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        r_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(),
                       weights=r_flux[ssort], cmap='gray',
                       range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])

        hdur = fits.PrimaryHDU(r[0] * weight_G)
        hdulr = fits.HDUList([hdur])
        file_name_r = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(
            theta_yz_1[0]) + '_rfilter.fits'
        hdulr.writeto('%s/%s/%s' % (str(save_directory),str(nout), file_name_r))

        hdur_non = fits.PrimaryHDU(r_non[0] * weight_G)
        hdulr_non = fits.HDUList([hdur_non])
        file_name_r_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(
            theta_yz_1[0]) + '_rfilter_non.fits'
        hdulr_non.writeto('%s/%s/%s' % (str(save_directory),str(nout), file_name_r_non))

    
        #np.save("./%s_%s_id" %(nout,idgal),star["id"][ssort])
        #np.save("./%s_%s_flux" %(nout,idgal),np.log10(r_flux[ssort]*1e3))
        #np.save("./%s_%s_flux_extinc" %(nout,idgal),np.log10(r_flux_extinc[ssort]*1e3))
        
        #save = open("%s_%s_flux.txt" %(nout,idgal),"a")
        #i = 0
        #while i < len(r_flux[ssort]):
        #    print(star["id"][ssort][i], np.log10(r_flux[ssort][i])+3, np.log10(r_flux_extinc[ssort][i])+3,file=save)
        #    i += 1
        #    try:
        #        print(np.round(100*i/len(r_flux[ssort])))
        #    except:
        #        pass
        #save.close()

    else:
        pass

    if saving == True and JHK == True:
        ### JHK ###
        sed.__init__(filter_system="JHK",info=gal.info)

        J_flux, J_flux_extinc = sed.get_flux(star, num=1, filter_name="J", cell=cell, colden=colden, extinction = True, info=gal.info)

        J = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=J_flux_extinc[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        J_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=J_flux[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])

        hduJ = fits.PrimaryHDU(J[0] * weight_B  )
        hdulJ = fits.HDUList([hduJ])
        file_name_J = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_Jfilter.fits'
        hdulJ.writeto('%s/%s' % (str(save_directory),file_name_J))

        hduJ_non = fits.PrimaryHDU(J_non[0] * weight_B)
        hdulJ_non = fits.HDUList([hduJ_non])
        file_name_J_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(
            theta_yz_1[0]) + '_Jfilter_non.fits'
        hdulJ_non.writeto('%s/%s' % (str(save_directory), file_name_J_non))



        H_flux, H_flux_extinc = sed.get_flux(star, num=1, filter_name="H", cell=cell, colden=colden, extinction = True, info=gal.info)

        H = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=H_flux_extinc[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        H_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=H_flux[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])

        hduH = fits.PrimaryHDU(H[0] * weight_G  )
        hdulH = fits.HDUList([hduH])
        file_name_H = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_Hfilter.fits'
        hdulH.writeto('%s/%s' % (str(save_directory),file_name_H))

        hduH_non = fits.PrimaryHDU(H_non[0] * weight_G)
        hdulH_non = fits.HDUList([hduH_non])
        file_name_H_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(
            theta_yz_1[0]) + '_Hfilter_non.fits'
        hdulH_non.writeto('%s/%s' % (str(save_directory), file_name_H_non))


        K_flux, K_flux_extinc = sed.get_flux(star, num=1, filter_name="K", cell=cell, colden=colden, extinction = True, info=gal.info)

        K = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=K_flux_extinc[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])
        K_non = plt.hist2d(star["y"][ssort], star["x"][ssort], bins=[npix, npix], norm=colors.LogNorm(), weights=K_flux[ssort], cmap='gray', range=[[drawlim_ybot,drawlim_yupp], [drawlim_xbot, drawlim_xupp]])

        hduK = fits.PrimaryHDU(K[0] * weight_R  )
        hdulK = fits.HDUList([hduK])
        file_name_K = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(theta_yz_1[0]) + '_Kfilter.fits'
        hdulK.writeto('%s/%s' % (str(save_directory),file_name_K))

        hduK_non = fits.PrimaryHDU(K_non[0] * weight_R)
        hdulK_non = fits.HDUList([hduK_non])
        file_name_K_non = str(nout) + '_' + str(idgal) + '_' + str(theta_xz_1[0]) + '_' + str(
            theta_yz_1[0]) + '_Kfilter_non.fits'
        hdulK_non.writeto('%s/%s' % (str(save_directory), file_name_K_non))


    
    else:
        pass



    ######################################################################
    print("save complete!")
    ############################################################3


    print("\n","x_real_span is ",x_real_span)
    print("\n","y_real_span is ",y_real_span)
    print("\n","z_real_span is ",z_real_span)

    print("\n","xspan is ",xspan*gal.info.boxtokpc)
    print("\n","yspan is ",yspan*gal.info.boxtokpc)
    print("\n","zspan is ",zspan*gal.info.boxtokpc)

    print("boxlim : histogram's range is",maxbin*gal.info.boxtokpc)


    ####################################################################################
    ######################## get additional properties #################################
    ####################################################################################

    ### GET SFR ###
    tc = utils.cosmology.Timeconvert(info=gal.info)
    starage = tc.time2gyr(star['time'],z_now = gal.info.zred)

    youngsort = np.where(
                            (star["x"] >= boxlim_xbot/gal.info.boxtokpc) & (star["x"] <= boxlim_xupp/gal.info.boxtokpc) &
                            (star["y"] >= boxlim_ybot/gal.info.boxtokpc) & (star["y"] <= boxlim_yupp/gal.info.boxtokpc) &
                            (starage > 0.0) & (starage <= 0.01)
                          )

    SFR = sum(star['m'][youngsort])



    ### GET COLD GAS MASS ###
    T_seperate = 6+0.25*np.log10(cell_rho)

    cold_gas_sort = np.where(  (np.log10(cell_T) < T_seperate)   )

    cold_gas_mass = cell_mgas[cold_gas_sort]
    print(np.log10(sum(cold_gas_mass)) )

    filename = "./properties_"+str(fixed_idgal)+"comp.txt"
    savefile = open(filename,"a")
    print(nout,idgal,sum(star['m']),sum(cell_mgas),sum(cold_gas_mass),SFR,VoS,DtT,
          np.sum(u_flux),np.sum(u_flux_extinc),          
          np.sum(g_flux),np.sum(g_flux_extinc),          
          np.sum(r_flux),np.sum(r_flux_extinc),          
          file=savefile)
    savefile.close()


    
    plt.close()

    return





##################################################

def auto_pixind(npix,real,pix_bot,pix_upp):
    center = npix/2
    stepsize = real/npix*2
    return (pix_bot-center)*stepsize, (pix_upp-center)*stepsize


##################################################

#nout_target,idgal_target = 599,2
fixlim = 20
pixlim = 1000

xbot = 0
xupp = pixlim
ybot = 0
yupp = pixlim
boxlim_xbot,boxlim_xupp = auto_pixind(pixlim,fixlim,xbot,xupp)
boxlim_ybot,boxlim_yupp = auto_pixind(pixlim,fixlim,ybot,yupp)

boxlim = fixlim
celllim = boxlim


draw_xbot = xbot
draw_xupp = xupp
draw_ybot = ybot
draw_yupp = yupp

drawlim_xbot = -fixlim
drawlim_xupp = +fixlim
drawlim_ybot = -fixlim
drawlim_yupp = +fixlim

drawpix = pixlim
drawlim_xbot,drawlim_xupp = auto_pixind(drawpix,fixlim,draw_xbot,draw_xupp)
drawlim_ybot,drawlim_yupp = auto_pixind(drawpix,fixlim,draw_ybot,draw_yupp)

#drawpix = 1000

#################################################

theta_xz_1 = [0]
theta_yz_1 = [0]

##################################################

nout_set = 100

dataroute = "/home/jangjk816/Project/Mock/Projenitors/"
dataname = str(nout_set)+"_projenitor.txt"
full_dataname = dataroute+dataname


dat = np.genfromtxt(full_dataname,dtype=[('nout','<i4'),('idgal','<i4')])
dat['nout'] += 39

k = 599
while k < 1000:
    
    try:
        load_gal(dat['nout'][599-k],dat['idgal'][599-k],
                theta_xz_1=theta_xz_1, theta_yz_1=theta_yz_1,
                      boxlim=boxlim, celllim=celllim,
                      boxlim_xbot=boxlim_xbot, boxlim_xupp=boxlim_xupp, boxlim_ybot=boxlim_ybot,
                      boxlim_yupp=boxlim_yupp,
                      drawlim_xbot=drawlim_xbot, drawlim_xupp=drawlim_xupp, drawlim_ybot=drawlim_ybot,
                      drawlim_yupp=drawlim_yupp, drawpix=drawpix,
                      saving=True, fuv=False, BVR=False, SDSS=True, JHK=False,
                      save_directory="./image/"+str(nout_set),fixed_idgal=nout_set)


    except:
        pass

    k += 1


#nout_set = 13
#fixed_idgal = 13
#theta_yz_1 = [0]

#load_gal(599,4,
#                theta_xz_1=theta_xz_1, theta_yz_1=theta_yz_1,
#                      boxlim=boxlim, celllim=celllim,
#                      boxlim_xbot=boxlim_xbot, boxlim_xupp=boxlim_xupp, boxlim_ybot=boxlim_ybot,
#                      boxlim_yupp=boxlim_yupp,
#                      drawlim_xbot=drawlim_xbot, drawlim_xupp=drawlim_xupp, drawlim_ybot=drawlim_ybot,
#                      drawlim_yupp=drawlim_yupp, drawpix=drawpix,
#                      saving=True, fuv=False, BVR=False, SDSS=True, JHK=False,
#                      save_directory="./image/"+str(nout_set),fixed_idgal=fixed_idgal)



###################################################

#76 ~ 523
#   0.8394253616634337
#   6.635994438423252 Gyr

#0 ~ 599
#   0.7021831598626218
#   7.356974742945828 Gyr


def continuous_figure(theta_yz_1_0):
    catalog = np.genfromtxt('13_progenitors.txt', dtype=[('snapnum','<f4'),('galnum','<f4'),('mass','<f4'),('idn','<f4')])
    l = len(catalog)-399
    while l < len(catalog) :
        try:
            load_gal(int(catalog['snapnum'][len(catalog)-l])+39,int(catalog['galnum'][len(catalog)-l]), theta_xz_1=theta_xz_1, theta_yz_1=theta_yz_1, boxlim=boxlim, celllim=celllim)
        except:
            pass
        #load_gal(int(catalog['snapnum'][len(catalog)-l])+39,int(catalog['galnum'][len(catalog)-l]), theta_xz_1=theta_xz_1, theta_yz_1=theta_yz_1, boxlim=boxlim, celllim=celllim)
        l = l +1



#continuous_figure(0)










