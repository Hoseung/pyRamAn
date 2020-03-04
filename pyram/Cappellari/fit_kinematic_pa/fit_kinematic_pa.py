"""
#############################################################################

 Copyright (C) 2005-2014, Michele Cappellari
 E-mail: cappellari_at_astro.ox.ac.uk

 Updated versions of the software are available from my web page
 http://purl.org/cappellari/software

 If you have found this software useful for your research,
 we would appreciate an acknowledgment to use of the
 `Method described in Appendix C of Krajnovic et al. (2006)'.

 This software is provided as is without any warranty whatsoever.
 Permission to use, for non-commercial purposes is granted.
 Permission to modify for personal or internal use is granted,
 provided this copyright and disclaimer are included unchanged
 at the beginning of the file. All other rights are reserved.

#############################################################################
+
 NAME:
   FIT_KINEMATIC_PA

 PURPOSE:
   Determine the global kinematic position angle of a
   galaxy with the method described in Appendix C of
   Krajnovic, Cappellari, de Zeeuw, & Copin 2006, MNRAS, 366, 787

 CALLING SEQUENCE:
   FIT_KINEMATIC_PA, Xbin, Ybin, Vel, AngleBest, AngleError, VelSyst, $
        /DEBUG, NSTEPS=nsteps, /QUIET

 INPUT PARAMETERS:
   XBIN, YBIN: vectors with the coordinates of the bins (or pixels)
       measured from the centre of rotation (typically the galaxy centre).
     - IMPORTANT: The routine will not give meaningful output unless 
       (X,Y)=(0,0) is an estimate of the centre of rotation.        
   VEL: measured velocity at the position (XBIN,YBIN). 
     - IMPORTANT: An estimate of the systemic velocity has to be already 
       subtracted from this velocity [e.g. VEL = VEL - median(VEL)]. 
       The routine will then provide in the output VELSYST a correction 
       to be added to the adopted systemic velocity.

 INPUT KEYWORDS:
   NSTEPS: number of steps along which the angle is sampled.
       Default is 361 steps which gives a 0.5 degr accuracy.
       Decrease this number to limit computation time during testing.

 OUTPUT PARAMETER:
   ANGLEBEST: kinematical PA. Note that this is the angle along which
       |Vel| is maximum (note modulus!). If one reverses the sense of
       rotation in a galaxy ANGLEBEST does not change. The sense of
       rotation can be trivially determined by looking at the map of Vel.
   ANGLEERROR: corresponding error to assign to ANGLEBEST.
   VELSYST: Best-fitting correction to the adopted systemic velocity 
       for the galaxy.
     - If the median was subtracted to the input velocity VEL before 
       the PA fit, then the corrected systemnic velocity will be 
       median(VEL)+VELSYST.

 REQUIRED ROUTINES:
   The following five additional routines are needed:
   - 1. CAP_SYMMETRIZE_VELFIELD and 2. CAP_RANGE: by Michele Cappellari
     (included in this FIT_KINEMATIC_PA distribution)
   - 3. SAURON_COLORMAP and 4. PLOT_VELFIELD: can be obtained from:
     http://purl.org/cappellari/idl#binning
   - 5. SIGRANGE: from IDL astro library http://idlastro.gsfc.nasa.gov/

 MODIFICATION HISTORY:
   V1.0 -- Created by Michele Cappellari, Leiden, 30 May 2005
   V1.1 -- Written documentation. MC, Oxford, 9 October 2007
   V1.11 -- Corrected handling of NSTEPS keyword. Thanks to Roland Jesseit.
       MC, Oxford, 17 October 2007
   V1.12 -- Force error to be smaller than 1/2 of the angular step size.
       MC, Oxford, 19 October 2007
   V1.13 -- Determine plotting ranges from velSym instead of vel.
       Thanks to Anne-Marie Weijmans. Leiden, 25 May 2008
   V1.14 -- Clarified that systemic velocity has to be subtracted from VEL.
       MC, Oxford, 31 March 2009
   V1.15 -- Overplot best PA on data. Some changes to the documentation.
       MC, Oxford, 14 October 2009
   V1.2 -- Includes error in chi^2 in the determination of angErr.
       Thanks to Davor Krajnovic for reporting problems.
       MC, Oxford, 23 March 2010
   V1.3 -- The program is two orders of magnitude faster, thanks to a 
       new cap_symmetrize_velfield routine. MC, Oxford, 8 October 2013
   V1.31 -- Uses CAP_RANGE routine to avoid potential naming conflicts.
       Uses TOLERANCE keyword of TRIANGULATE to try to avoid IDL error
       "TRIANGULATE: Points are co-linear, no solution." 
       MC, Oxford, 2 December 2013
   V2.0 -- Translated from IDL into Python. MC, Oxford, 10 April 2014
   V2.01 -- Support both Python 2.6/2.7 and Python 3.x. MC, Oxford, 25 May 2014

"""

from __future__ import print_function

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from cap_symmetrize_velfield import symmetrize_velfield
from cap_plot_velfield import plot_velfield

##############################################################################

def fit_kinematic_pa(x, y, vel, debug=False, nsteps=361, 
                     quiet=False, plot=True, dvel=None):

    if dvel is None:
        dvel = vel*0 + 10.0 # Adopt here constant 10 km/s errors!
    
    nbins = x.size
    n = nsteps
    angles = np.linspace(0, 180, n) # 0.5 degrees steps by default
    chi2 = np.empty_like(angles)
    for j, ang in enumerate(angles):
        velSym = symmetrize_velfield(x, y, vel, sym=1, pa=ang)
        chi2[j] = np.sum(((vel-velSym)/dvel)**2)
        if debug:
            print('Ang, chi2/DOF:', ang, chi2[j]/nbins)
            plot_velfield(x, y, velSym)
            plt.pause(0.01)
    k = np.argmin(chi2)
    angBest = angles[k]
    
    # Compute fit at the best position
    #
    velSym = symmetrize_velfield(x, y, vel, sym=1, pa=angBest)
    if angBest < 0:
        angBest += 180
    
    # 3sigma confidence limit, including error on chi^2
    #
    f = chi2 - chi2[k] <= 9 + 3*np.sqrt(2*nbins)
    if f.sum():
        angErr = (np.max(angles[f]) - np.min(angles[f]))/2.0
        if angErr >= 45:
            good = np.degrees(np.arctan(np.tan(np.radians(angles[f]))))
            angErr = (np.max(good) - np.min(good))/2.0
    else:
        angErr = max(0.5, (angles[1]-angles[0])/2.0)
    
    angErr = angErr.clip(max(0.5, (angles[1]-angles[0])/2.0)) # Force errors to be larger than 0.5 deg
    vSyst = np.median(vel - velSym)
    
    if not quiet:
        print('  Kin PA:', angBest, ' +/- ', angErr, ' (3*sigma error)')
        print('Velocity Offset:', vSyst)
    
    # Plot results
    #
    if plot:    
    
        mn, mx = stats.scoreatpercentile(velSym, [2.5, 97.5])
        mx = min(mx, -mn)
        plt.subplot(121)
        plot_velfield(x, y, velSym, vmin=-mx, vmax=mx) 
        plt.title('Symmetrized')
        
        plt.subplot(122)
        plot_velfield(x, y, vel-vSyst, vmin=-mx, vmax=mx) 
        plt.title('Data and best PA')
        rad = np.sqrt(np.max(x**2 + y**2))
        ang = [0,np.pi] + np.radians(angBest)
        plt.plot(rad*np.cos(ang), rad*np.sin(ang), '--', linewidth=3) # Zero-velocity line
        plt.plot(-rad*np.sin(ang), rad*np.cos(ang), linewidth=3) # Major axis PA

    return angBest, angErr, vSyst
    
##############################################################################

# Usage example for fit_kinematic_pa()

if __name__ == '__main__':
    
    xbin, ybin = np.random.uniform(low=[-30,-20], high=[30,20], size=(100,2)).T
    inc = 60.                       # assumed galaxy inclination
    r = np.sqrt(xbin**2 + (ybin/np.cos(np.radians(inc)))**2) # Radius in the plane of the disk
    a = 40                          # Scale length in arcsec
    vr = 2000*np.sqrt(r)/(r+a)      # Assumed velocity profile
    vel = vr * np.sin(np.radians(inc))*xbin/r # Projected velocity field

    plt.clf()
    fit_kinematic_pa(xbin, ybin, vel)
    