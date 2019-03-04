from collections import namedtuple

Cosmo = namedtuple("cosmo_params", ["Om", "Ol", "Ob", "H0", "sig8", "nspec"])
Planck15 = Cosmo(0.3089, 0.6911, 0.0460, 67.74, 0.8159, 0.9667)
Planck18 = Cosmo(0.3111, 0.6889, 0.0488, 67.66, 0.8102, 0.9655)

Consts = namedtuple("Constants", ["G", "Msun_in_g", "pc_in_cm"])
# All in cgs
consts = Consts( 6.67408e-8, 1.989e33, 3.08567758e18)
