import ramski 
import numpy as np
import pyram
from pyram.utils.cosmology import Timeconvert as TC


YZiCS = False
NH = True

if YZiCS:
    centers = (0.51237148,  0.68784678,  0.45724899)# JP's pick 
    nout = 187 
    ramski.ramski_v4("/home/hoseung/Work/data/YZiCS/29172/snapshots/", 
                 "./output/", 
                 "/home/hoseung/Work/data/YZiCS/29172/halo/DM/tree.dat", 
                 centers, nout)

if NH:
    wdir = '/home/hopung/Work/NH/'
    outdir = './'
    tdir = '/storage1/NewHorizon/GalaxyMaker/gal/tree.dat',
    centers = (0.49977215140452347, 0.4922276662830468, 0.50730382849124867)
    nout = 509
    ramski.ramski_v4(wdir+'snapshots/',
                    outdir,
                    './',
                    centers,nout)

## Particle
s=pyram.load.sim.Sim(nout, base=wdir)
radius = 0.00011

s.set_ranges([[centers[0]-radius, centers[0]+radius],
              [centers[1]-radius, centers[1]+radius],
              [centers[2]-radius, centers[2]+radius]])

# amr2cell only reads reaf cells. So, ref == 0 for all cell.
#s.add_hydro(ref=False, additional_field=[{"name":"ind", "dtype":"<i4", "dim":1}])

s.add_part(ptypes=["star id pos mass vel metal age"])
star = s.part.star
tc = TC(s.info)
age = tc.time2gyr(times = star["time"], z_now = s.info.zred)

# write part.ski
"""
http://www.skirt.ugent.be/skirt9/class_text_in_file.html
There must be a header line for each column with the following contents from left to right:

the hash character indicating a header line
the word "column" (case insensitive)
the one-based column number (optional)
a colon
a description that does not contain parenthesis (optional)
a unit string between parenthesis (may be "(1)" or "()" for dimensionless quantities)
"""

# Star particle for testing SED family
# Column 1: position x (pc)
# Column 2: position y (pc)
# Column 3: position z (pc)
# Column 4: smoothing length (pc)
# Column 5: mass (Msun)
# Column 6: metallicity (1)
# Column 7: age (Gyr)
# 0 0 0 1   10  0.04  5

pos_pc = (star["pos"] - centers) *s.info.boxtokpc * 1e3

smooth = 25 # arbitrary

with open("part_ramski.txt", "w") as fpart:
    fpart.write("""# Column 1: position x (pc)
# Column 2: position y (pc)
# Column 3: position z (pc)
# Column 4: smoothing length (pc)
# Column 5: mass (Msun)
# Column 6: metallicity (1)
# Column 7: age (Gyr)
""")
    for (pos, m, z, a) in zip(pos_pc, star["m"] * s.info.msun, star["metal"], age):
        fpart.write("{} {} {} {} {} {} {}\n".format(pos[0], pos[1], pos[2], smooth, m, z, a))
       
print("done") 
