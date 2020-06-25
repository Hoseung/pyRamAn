import numpy as np
import pickle
import pyram

from pyram.tree import tmtree as tmt
from pyram.utils.cosmology import Timeconvert as TC
from pyram.mock import ramski # requires ramski_module.so
from pyram import galaxymodule as gmo
import pts.simulation as sm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def write_ski(repo, gid, nout, skifile, arr):
    inclination, azimuth, roll, fovx, fovy, minX, maxX, minY, maxY, minZ, maxZ = arr

    # load values from JP's value

    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'minX', minX)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'maxX', maxX)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'minY', minY)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'maxY', maxY)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'minZ', minZ)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'maxZ', maxZ)

    # PTS doesn't support modifying values if attribute has multiple elments.
    elems = skifile._tree.xpath('//MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/FullInstrument')

    inst_fc = elems[0] # Instrument "faceon"
    inst_fc.set('instrumentName', 'faceon')
    inst_fc.set('inclination', f'{inclination:.5f} deg')
    inst_fc.set('azimuth', f'{azimuth:.5f} deg')
    inst_fc.set('roll', f'{roll:.5f} deg')

    #edge-on
    if inclination < 90: 
        #print("Warning... inclination =", inclination)
        inclination0 = inclination + 90
    else:
        inclination0 = inclination - 90
        #azimuth
    for i, inst_ed in enumerate(elems[1:]):
        inclination = inclination0 - 9 + 3*i 
        if inclination < 0 or inclination > 180:
            continue
        inst_ed.set('instrumentName', f'edgeon{i}')

        inst_ed.set('inclination', f'{inclination:.5f}')
        inst_ed.set('azimuth', f'{azimuth:.5f}')
        inst_ed.set('roll', f'{roll:.5f}')

        inst_ed.set('fieldOfViewX', f'{fovx} pc')
        inst_ed.set('fieldOfViewY', f'{fovy} pc')

    skifile.saveTo(repo+f"g{gid}_{nout}.ski")
    print(f"Done {gid}")




def make_ram_input(wdir, nout, gal, dir_out='./', fov=40, smooth = 25, nvec_frac=-1, fsave_angle=None, plot_stellar=False):
    """
    fov in kpc. 

    """
    centers, radius, gid=  (gal['x'],gal['y'],gal['z']), gal['r'], gal['id']
    
    arr=ramski.ramski_v4(wdir+'snapshots/',
                    dir_out,
                    centers,
                    nout,
                    radius,
                    gid)

    incl, azim, roll, fovx, fovy, minx, maxx, miny, maxy, minz, maxz = arr
    # Note that I will not use incl, azim from this measurement.
    print("min, max", minx, maxx, miny, maxy, minz, maxz)
    print("radius", radius)

    ## Particle
    s=pyram.load.sim.Sim(nout, base=wdir)

    s.set_ranges([[centers[0]-radius, centers[0]+radius],
                  [centers[1]-radius, centers[1]+radius],
                  [centers[2]-radius, centers[2]+radius]])

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

    ### Write particles
    with open(dir_out+ f"{gid:05d}/{nout:05d}/part_ramski.txt", "w") as fpart:
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

    if nvec_frac > 0:
        gal = gmo.galaxy.Galaxy()
        gal.info = s.info
        gal.star = s.part.star

        vcen = np.mean(gal.star['vel'], axis=0)
        gal.star['vel'] -= vcen

        print("# {} stars in total".format(len(gal.star)))

        dist = np.sqrt(np.einsum("...i,...i", gal.star["pos"],gal.star["pos"]))
        ind_close = np.argsort(dist)[:int(nvec_frac*len(gal.star))]
        #close_stars = gal.star[ind_close]

        params_ang = np.mean(gal.star[ind_close]['m'] * np.cross(gal.star[ind_close]["pos"], gal.star[ind_close]["vel"]).T, axis=1)
        #close_stars = None
        print("params_ang", params_ang)

        lvec = np.sum(np.cross(gal.star['pos'][ind_close],
                               gal.star['vel'][ind_close]), axis=0)
        nvec = lvec/np.linalg.norm(lvec)

        theta = np.degrees(np.arccos(nvec[2]))
        if nvec[1] > 0:
            phi = np.degrees(np.arccos(nvec[0]/np.sqrt(nvec[0]**2 + nvec[1]**2)))
        else:
            phi = np.degrees(-np.arccos(nvec[0]/np.sqrt(nvec[0]**2 + nvec[1]**2)))
        if fsave_angle is not None: fsave_angle.write(f"{nout} {gid} {theta:.3f}  phi={phi:.3f} \n")
        arr[:3] = theta, phi, 0
        print(f"inclination after {theta:.3f}")
        print(f"azimuth after {phi:.3f}")

        if plot_stellar:
            fig, axs = plt.subplots(2,2)
            fig.set_size_inches(8,8)
            axs[0,0].hist2d(gal.star["x"], gal.star["y"], bins=400, norm=LogNorm())
            gal.meta.nvec = nvec
            gal.reorient(dest=[0,0,1]) 
            axs[0,1].hist2d(gal.star["x"], gal.star["y"], bins=400, norm=LogNorm())
            axs[1,0].hist2d(gal.star["x"], gal.star["z"], bins=400, norm=LogNorm())
            for ax in axs.ravel():
                ax.set_aspect("equal")
            plt.savefig(dir_out+ f"{gid:05d}_stellar_maps.png", dpi=200)
    
    # update .ski immediately
    fn_template = './template_zubko.ski'
    skifile = sm.SkiFile(fn_template)
    repo = dir_out+ f"{gid:05d}/{nout:05d}/"#faceon_redshift_"
    write_ski(repo, gid, nout, skifile, arr)


    print("done") 


if __name__ == "__main__":
    nout=906
    gal_cat = np.array([13, 6.854e+10, 0.48733889, 0.47836283, 0.49842872, 1.5930e-04])
    #gcats = np.genfromtxt(f"./centers_{nout}.txt",
    #                       dtype=[('id','int'), ('x','float'),('y','float'),('z','float')])

    radius = 20/(100*1e3/0.704)
    with open(f"rot_angles_{nout}.txt", "w") as fsave:
        for gal in gcats:
            print(gal)
            make_ram_input("./", nout, gal,
                        dir_out=f'./gals_{nout}/', nvec_frac=0.4, fsave_angle=fsave)


