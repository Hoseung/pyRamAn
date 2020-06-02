import numpy as np
import pickle
import pyram

from pyram.tree import tmtree as tmt
from pyram.utils.cosmology import Timeconvert as TC
from pyram.mock import ramski # requires ramski_module.so
from pyram import galaxymodule as gmo

def make_ram_input(wdir, nout, centers, radius, gid, fov=40, smooth = 25):
    ramski.ramski_v4(wdir+'snapshots/',
                    './',
                    centers,
                    nout,
                    fov,
                    gid)

    ## Particle
    s=pyram.load.sim.Sim(nout, base=wdir)

    s.set_ranges([[centers[0]-radius, centers[0]+radius],
                  [centers[1]-radius, centers[1]+radius],
                  [centers[2]-radius, centers[2]+radius]])

    # amr2cell only reads reaf cells. So, ref == 0 for all cell.
    #s.add_hydro(ref=False, additional_field=[{"name":"ind", "dtype":"<i4", "dim":1}])

    s.add_part(ptypes=["star id pos mass vel metal age"])
    star = s.part.star
    tc = TC(s.info)
    age = tc.time2gyr(times = star["time"], zred_now = s.info.zred)
    age = -1*age*1e9
    # Make gal out of part 

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
    pos_pc = (star["pos"] - centers) *s.info.boxtokpc * 1e3

    with open(f"./{gid:05d}/{nout:05d}/part_ramski.txt", "w") as fpart:
        fpart.write("""# Column 1: position x (pc)
    # Column 2: position y (pc)
    # Column 3: position z (pc)
    # Column 4: smoothing length (pc)
    # Column 5: mass (Msun)
    # Column 6: metallicity (1)
    # Column 7: age (Gyr)\n
    """)
        for (pos, m, z, a) in zip(pos_pc, star["m"] * s.info.msun, star["metal"], age):
            fpart.write("{} {} {} {} {} {} {}\n".format(pos[0], pos[1], pos[2], smooth, m, z, a))

    print("done") 


nout=874
#gal_cat = np.array([13, 6.854e+10, 0.48733889, 0.47836283, 0.49842872, 1.5930e-04])
gcats = np.genfromtxt(f"./centers_{nout}.txt",
                       dtype=[('id','int'), ('x','float'),('y','float'),('z','float')])

gal = gcats[50]
radius = 20/(100*1e3/0.704)
make_ram_input("./", nout, (gal['x'],gal['y'],gal['z']), radius, gal['id'])


