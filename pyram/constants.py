from collections import namedtuple

"""
Provide as complete as possible list of entries as namedtuple is immutable after created.
"""

Cosmo = namedtuple("cosmo_params", ["Om", "Ol", "Ob", "H0", "sig8", "nspec"])
Planck15 = Cosmo(0.3089, 0.6911, 0.0460, 67.74, 0.8159, 0.9667)
Planck18 = Cosmo(0.3111, 0.6889, 0.0488, 67.66, 0.8102, 0.9655)

Consts = namedtuple("Constants", ["G", "G_mks", "Msun_in_g", "pc_in_cm", "Mpc_in_km"])
# All in cgs
consts = Consts( 6.67408e-8, 6.67408e-11, 1.989e33, 3.08567758e18, 3.08567758e19)

# some units in cgs
Msol = 2E33 # Msol = 1.9885E33

# parsecs: ramses unit
Mpc = 3.08E24
kpc = 3.08E21
pc = 3.08E18
ly = 9.4607E17
km = 1E5

Gyr = 3.154E16
Myr = 3.154E13
yr = 3.154E7

# some constants in cgs
k_B = 1.38064852E-16 # Boltzmann constant
m_H = 1.6737236E-24 # hydrogen atomic mass

# others
gamma = 1.6666667


# in cgs unit
kpc = 3.08e21
twopi = 6.2831853e0
hplanck = 6.6262000e-27
eV = 1.6022000e-12
kB = 1.38e-16
clight = 2.9979250e+10
Gyr = 3.1536000e+16
X = 0.76
Y = 0.24
rhoc = 1.8800000e-29
mH = 1.6600000e-24
mu_mol = 1.2195e0
G = 6.67e-8
m_sun = 1.98892e33