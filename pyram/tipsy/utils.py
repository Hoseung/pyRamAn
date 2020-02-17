#Om = 0.3089
from astropy.cosmology import Planck15
import astropy.units as u

pp = Planck15

rho_crit = pp.critical_density0

## How do I calculate rho_crit?? 
# Look at astropy/cosmology/core.py

def get_Msol(rho_crit, lbox_in_kpc = 9e4):
    Msun = 1.989e33
    kpc_in_cm = 3.08567758e21
    # 147798.126
    V_in_cc = (kpc_in_cm*lbox_in_kpc)**3
    
    dMsolUnit =  rho_crit / (Msun / V_in_cc)
    return dMsolUnit
