#!/usr/bin/env python

"""
    This example illustrates how to obtain an MGE fit from a galaxy image
    using the mge_fit_sectors package and how to verify the results.

    V1.0.0: Translated NGC4342 example from the corresponding IDL version.
        Michele Cappellari, Aspen Airport, 8 February 2014
    V1.0.1: Fixed incorrect offset in high-res contour plot.
        Use arcsec pixel scale. MC, Oxford, 18 September 2014
    V1.1.0: Translated M32 example of fitting to two images simultaneously.
        MC, Oxford, 18 June 2015
    V1.1.1: Support both Pyfits and Astropy to read FITS files.
        MC, Oxford, 23 October 2015

"""

import numpy as np
import matplotlib.pyplot as plt
try:
    import pyfits
except:
    from astropy.io import fits as pyfits

from sectors_photometry import sectors_photometry
from mge_print_contours import mge_print_contours
from find_galaxy import find_galaxy
from mge_fit_sectors import mge_fit_sectors

#----------------------------------------------------------------------------

def fit_m32():
    """
    This procedure reproduces Figures 6-7 in Cappellari (2002).

    This example illustrates how to fit multiple images together.
    We model an HST/WFPC2/F814W image of M32 and an I-band ground-based one.

    """

    scale1 = 0.0455     # (arcsec/pixel) PC1. This is used as scale and flux reference!
    scale2 = 0.549      # (arcsec/pixel) ING (Peletier 1993)
    scaleRatio = scale1/scale2
    fluxRatio = 0.9579  # = flux1/flux2 (ratio of observed counts/pixel at give radius)

    # IMPORTANT: use the *same* eps to run sectors_photometry on both images
    eps = 0.2

    # Perform photometry on the HST/WFPC2/F814W (I-band) image
    # The geometric parameters below were obtained using my FIND_GALAXY program
    ang1 = 165.4
    xc1 = 377
    yc1 = 314

    file = "m32_f814w_pc.fits"
    hdu = pyfits.open(file)
    img1 = hdu[0].data
    img1 -= 4.48  # subtract sky

    plt.clf()
    s1 = sectors_photometry(img1, eps, ang1, xc1, yc1, minlevel=0, plot=1)
    plt.pause(0.01)  # Allow plot to appear on the screen

    # Perform photometry on Peletier (1993) ING/I-band image
    # The geometric parameters below were obtained using my FIND_GALAXY program
    ang2 = 22.9
    xc2 = 376
    yc2 = 184

    file = "m32_i.fits"
    hdu = pyfits.open(file)
    img2 = hdu[0].data
    img2 -= 32.0  # subtract sky, determined to make outer profile asymptotically power-law

    plt.clf()
    s2 = sectors_photometry(img2, eps, ang2, xc2, yc2, minlevel=5.0, plot=1)
    plt.pause(0.01)  # Allow plot to appear on the screen
    s2.radius /= scaleRatio  # Bring all radii and fluxes on the same scale
    s2.counts *= fluxRatio

    # Exclude pixels at small radii (<3") in Peletier's image to avoid
    # PSF effects and merges the profiles of the different images.
    # The HST image is used as flux and spatial scale reference,
    # the ground-based data were simply scaled to the HST units.
    w = s2.radius > 3/scale1
    radius = np.append(s1.radius, s2.radius[w])
    angle = np.append(s1.angle, s2.angle[w])
    counts = np.append(s1.counts, s2.counts[w])

    # The PSF needs to be the one for the high-resolution image used in the centre.
    # Here this is the WFPC2/PC1/F814W image (we use a Gaussian PSF for simplicity)
    sigmaPSF = 0.8
    ngauss = 11

    # Do the actual MGE fit
    m = mge_fit_sectors(radius, angle, counts, eps, ngauss=ngauss, sigmaPSF=sigmaPSF,
                        scale=scale1, plot=1, linear=0, qbounds=[0.3, 0.85])
    plt.pause(0.01)  # Allow plot to appear on the screen

    plt.subplot(121)
    # Plot MGE contours on top of the HST image
    mge_print_contours(img1, ang1, xc1, yc1, m.sol, scale=scale1,
                       binning=4, sigmapsf=sigmaPSF, magrange=9)

    # Scale the solution parameters to the ground-based image
    m.sol[0, :] *= scaleRatio**2/fluxRatio  # Gaussians total counts
    m.sol[1, :] *= scaleRatio               # sigma of the Gaussians
    sigmaPSF = 1/2.35  # seeing FWHM = 1.0 arcsec for the ground based image

    plt.subplot(122)
    # Plot MGE contours on top of the ground-based image
    mge_print_contours(img2, ang2, xc2, yc2, m.sol, scale=scale2,
                       binning=4, sigmapsf=sigmaPSF, magrange=9.5)
    plt.pause(0.01)  # Allow plot to appear on the screen

#----------------------------------------------------------------------------

def fit_ngc4342():
    """
    This procedure reproduces Figures 8-9 in Cappellari (2002)
    This example illustrates a simple MGE fit to one single HST/WFPC2 image.

    """

    file = "ngc4342_f814w_pc.fits"

    hdu = pyfits.open(file)
    img = hdu[0].data

    skylev = 0.55   # counts/pixel
    img -= skylev   # subtract sky
    scale = 0.0455  # arcsec/pixel
    minlevel = 0.2  # counts/pixel
    ngauss = 12

    # Here we use an accurate four gaussians MGE PSF for
    # the HST/WFPC2/F814W filter, taken from Table 3 of
    # Cappellari et al. (2002, ApJ, 578, 787)

    sigmaPSF = [0.494, 1.44, 4.71, 13.4]      # In PC1 pixels
    normPSF = [0.294, 0.559, 0.0813, 0.0657]  # total(normPSF)=1

    # Here we use FIND_GALAXY directly inside the procedure. Usually you may want
    # to experiment with different values of the FRACTION keyword, before adopting
    # given values of Eps, Ang, Xc, Yc.
    plt.clf()
    f = find_galaxy(img, fraction=0.04, plot=1)
    plt.pause(0.01)  # Allow plot to appear on the screen

    # Perform galaxy photometry
    plt.clf()
    s = sectors_photometry(img, f.eps, f.theta, f.xpeak, f.ypeak,
                           minlevel=minlevel, plot=1)
    plt.pause(0.01)  # Allow plot to appear on the screen

    # Do the actual MGE fit
    m = mge_fit_sectors(s.radius, s.angle, s.counts, f.eps,
                        ngauss=ngauss, sigmaPSF=sigmaPSF, normPSF=normPSF,
                        scale=scale, plot=0, bulge_disk=0, linear=0)

    # Show contour plots of the results
    plt.subplot(121)
    mge_print_contours(img, f.theta, f.xpeak, f.ypeak, m.sol, scale=scale,
                       binning=7, sigmapsf=sigmaPSF, normpsf=normPSF, magrange=9)

    # Extract the central part of the image to plot at high resolution.
    # The MGE is centered to fractional pixel accuracy to ease visual comparson.

    n = 50
    img = img[f.xpeak-n:f.xpeak+n, f.ypeak-n:f.ypeak+n]
    xc, yc = n - f.xpeak + f.xmed, n - f.ypeak + f.ymed
    plt.subplot(122)
    mge_print_contours(img, f.theta, xc, yc, m.sol,
                       sigmapsf=sigmaPSF, normpsf=normPSF, scale=scale)
    plt.pause(0.01)  # Allow plot to appear on the screen

#----------------------------------------------------------------------------

if __name__ == '__main__':

    print("\nFitting M32---------------------------------------------\n")
    fit_m32()

    print("\nFitting NGC4342-----------------------------------------\n")
    fit_ngc4342()
