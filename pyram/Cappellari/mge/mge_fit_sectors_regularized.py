################################################################################
#
# Copyright (C) 2013-2015, Michele Cappellari
# E-mail: michele.cappellari_at_physics.ox.ac.uk
#
# Updated versions of the software are available from my web page
# http://purl.org/cappellari/software
#
# This software is provided as is without any warranty whatsoever.
# Permission to use, for non-commercial purposes is granted.
# Permission to modify for personal or internal use is granted,
# provided this copyright and disclaimer are included unchanged
# at the beginning of the file. All other rights are reserved.
#
################################################################################
#
# This program is a wrapper for mge_fit_sectors procedure and it accepts all
# keyword of that program. One should look at the documentation of mge_fit_sectors
# for usage details.
#
# The wrapper implements the method described in Section 2.2.2 of
# Cappellari (2002, MNRAS, 333, 400) to "regularize" an MGE model by restricting
# the allowed range in qObs of the Gaussians until the fit becomes unacceptable.
# In this way one ensures that the permitted galaxy inclinations are not being
# artificially restricted to a smaller range than allowed by the data.
#
# The detailed approach implemented here is the one used for the MGE fits of
# the galaxies in the Atlas3D project and described in Section 3.2 of
# Scott et al. (2013, MNRAS, 432, 1894).
#
# The intended usage of this wrapper is the following:
#   1. First perform a standard MGE fit with mge_fit_sectors;
#   2. Once all parameters (e.g. PA, eps, centre, sky subtraction) are OK and the
#      fit looks good, simply rename "mge_fit_sectors" in your script into
#      "mge_fit_sectors_regularized" (and import the module) to cleanup the final solution.
#      This is because this wrapper calls mge_fit_sectors repeatedly, taking much longer,
#      so it is not useful to run it until all input parameters are settled.
#
# VERSION HISTORY:
#   V1.0.0: Michele Cappellari, Oxford, 22 January 2013
#   V1.0.1: Fixed program stop when (qmin==qmax). Thanks to Silvia Posacki (Bologna)
#       for reporting the problem and the solution. MC, Oxford, 9 May 2013
#   V2.0.0: Converted from IDL into Python. MC, Oxford, 27 March 2015
#   V2.0.1; Removed truncation of input eps in mge_fit_sectors.
#       MC, Atlantic Ocean, 28 March 2015
#   V2.0.2: Cleaned up loop. MC, Oxford, 30 May 2015
#
################################################################################

import numpy as np

from Cappellari.mge.mge_fit_sectors import mge_fit_sectors

#----------------------------------------------------------------------------

class mge_fit_sectors_regularized(object):

    def __init__(self, radius, angle, counts, eps, qbounds=[0, 1], **kwargs):

        qmin, qmax = qbounds
        nq = int(np.ceil((qmax - qmin)/0.05) + 1)  # Adopt step <= 0.05 in qObs
        qrange = np.linspace(qmin, qmax, nq)
        bestnorm = np.inf
        frac = 1.1  # Allowed fractional increase in ABSDEV

        for j in range(nq - 1):
            qmin = qrange[j]
            m = mge_fit_sectors(radius, angle, counts, eps, qbounds=[qmin, qmax], **kwargs)
            absdev = m.absdev
            print('(minloop) qbounds=%6.4f %6.4f' % (qmin, qmax))
            if absdev > bestnorm*frac:
                jbest = j - 1
                qmin = qrange[jbest]
                break  # stops if error increases more than frac
            else:
                jbest = j
                bestnorm = min(bestnorm, absdev)
                self.sol = m.sol

        for k in range(nq - 2, jbest, -1):
            qmax = qrange[k]
            m = mge_fit_sectors(radius, angle, counts, eps, qbounds=[qmin, qmax], **kwargs)
            absdev = m.absdev
            print('(maxloop) qbounds=%6.4f %6.4f' % (qmin, qmax))
            if absdev > bestnorm*frac:
                qmax = qrange[k + 1]
                break  # stops if error increases more than frac
            else:
                bestnorm = min(bestnorm, absdev)
                self.sol = m.sol

        print('Final qbounds=%6.4f %6.4f' % (qmin, qmax))

#----------------------------------------------------------------------------
