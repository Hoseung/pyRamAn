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

def write_ski(repo, gid, nout, skifile, arr, params, pixel_scale=35):

    inclination, azimuth, roll, fovx, fovy, minX, maxX, minY, maxY, minZ, maxZ = arr
    if np.sum(arr[:3]) > 0:
    # so that rotated stars all in the FoV. 
        fovx *= 1.7
        fovy *= 1.7 
    else:
        fovx*=2
        fovy*=2
        print("No need to increase FoV")
    # load values from JP's value
    skifile.setIntAttribute('//MonteCarloSimulation', 'numPackets', params[0]['nphoton'])
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'minX', minX)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'maxX', maxX)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'minY', minY)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'maxY', maxY)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'minZ', minZ)
    skifile.setFloatAttribute('//MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium', 'maxZ', maxZ)

    # PTS doesn't support modifying values if attribute has multiple elments.
    elems = skifile._tree.xpath('//MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/FullInstrument')

    #print("length of params,", len(params))
    for param, inst in zip(params, elems[:len(params)]):
        inst.set('instrumentName', param['name'])
        inst.set('inclination', '{:.5f} deg'.format(inclination + param['incli_off']))
        inst.set('azimuth', '{:.5f} deg'.format(azimuth + param['azim_off']))
        inst.set('roll', '{:.5f} deg'.format(roll + param['roll_off']))
    
        inst.set('fieldOfViewX', f'{fovx} pc')
        inst.set('fieldOfViewY', f'{fovy} pc')

        npixx = npixy = int(fovx / param['pixel_scale'])
        inst.set('numPixelsX', f'{npixx}')
        inst.set('numPixelsY', f'{npixy}')

    skifile.saveTo(repo+f"g{gid}_{nout}.ski")
    print(f"Done {gid}")


def make_ram_input(wdir, nout, gcat, fsave,
                          dir_out='./',
                          ski_params=None,
                          smooth = 50,
                          nvec_frac=-1,
                          plot_stellar=False,
                          more_props=False,
                          r_fixed=None):
    """
    fov in kpc. 

    """
    if ski_params == None:
        ski_params=[dict(name='face_on',  pixel_scale=50, nphoton=5e7, incli_off = 0),
                dict(name='edge_on',  pixel_scale=50, nphoton=5e7, incli_off = 90)]

    centers, radius, gid = (gcat['x'],gcat['y'],gcat['z']), gcat['r'], gcat['id']

    s=pyram.load.sim.Sim(nout, base=wdir)
    radius_in_kpc = radius*s.info.boxtokpc
    if r_fixed is not None:
        radius_in_kpc = r_fixed

    arr=ramski.ramski_v4(wdir+'snapshots/',
                    dir_out,
                    centers,
                    nout,
                    radius_in_kpc,
                    gid)
    # fov in kpc -> in pc
    arr[3:5] *= 1e3
    incl, azim, roll, fovx, fovy, minx, maxx, miny, maxy, minz, maxz = arr
    # Note that I will not use incl, azim from this measurement.
    print("min, max", minx, maxx, miny, maxy, minz, maxz)
    print("radius", radius_in_kpc)
    print("FoV", fovx, fovy)
    print(minx, maxx, miny, maxy, minz, maxz)

    print("Xrange", centers[0]-radius, centers[0]+radius)
    ## Particle
    s.set_ranges([[centers[0]-radius, centers[0]+radius],
                  [centers[1]-radius, centers[1]+radius],
                  [centers[2]-radius, centers[2]+radius]])

    s.add_part(ptypes=["star id pos mass vel metal age"])
    star = s.part.star
    star['m'] *= s.info.msun
    tc = TC(s.info)
    age = tc.time2gyr(times = star["time"], zred_now = s.info.zred)
    age = s.info.tGyr - age
    age[age<0] = 1e-7
    age *= 1e9
    print("Age min, max", min(age), max(age))
    # Make gal out of part 

    # write part.ski
    """
    http://www.skirt.ugent.be/skirt9/class_text_in_file.html
    There must be a header line for each column with the following contents from left to right:
    hash character indicating a header line
    the word "column" (case insensitive)
    the one-based column number (optional)
    a colon
    a description that does not contain parenthesis (optional)
    a unit string between parenthesis (may be "(1)" or "()" for dimensionless quantities)
    """
    # Why do I need a new center when GalaxyMaker's center is precise enough? 
    #centers = np.mean(star['pos'], axis=0)
    #print("new center", centers)
    pos_pc = (star["pos"] - centers) *s.info.boxtokpc * 1e3

    ### Write particles
    with open(dir_out+ f"{gid:05d}/{nout:05d}/part_ramski.txt", "w") as fpart:
        fpart.write("""# Column 1: position x (pc)
    # Column 2: position y (pc)
    # Column 3: position z (pc)
    # Column 4: smoothing length (pc)
    # Column 5: mass (Msun)
    # Column 6: metallicity (1)
    # Column 7: age (yr)\n
    """)
        for (pos, m, z, a) in zip(pos_pc, star["m"], star["metal"], age):
            fpart.write("{} {} {} {} {} {} {}\n".format(pos[0], pos[1], pos[2], smooth, m, z, a))

    gal = gmo.galaxy.Galaxy()
    gal.info = s.info
    gal.star = s.part.star
    vcen = np.mean(gal.star['vel'], axis=0)
    gal.star['vel'] -= vcen
    gal.star['pos'] -=centers
    gal.star['pos'] *= gal.info.boxtokpc

    print("# {} stars in total".format(len(gal.star)))
    print("pos min max", np.min(gal.star['pos'], axis=0), np.max(gal.star['pos'],axis=0))

    if nvec_frac > 0:
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
        alpha = 0

    else:
        nvec = np.array([0,0,1])
        theta, phi, alpha = 0,0,0

    if not more_props:
        fsave.write(f"{nout} {gid} {theta:.6f}  {phi:.6f} {centers[0]:.6f} {centers[1]:.6f} {centers[2]:.6f} {nvec[0]:.6f} {nvec[1]:.6f} {nvec[2]:.6f} \n")
    else:
        new_star = np.sum(gal.star['m'][age < 1e8])
        fsave.write(f"{nout} {gid} {theta:.4f}  {phi:.4f} {centers[0]:.5f} {centers[1]:.5f} {centers[2]:.5f} {nvec[0]:.5f} {nvec[1]:.5f} {nvec[2]:.5f} " + 
                    "{:.5f} {:.5f} {:.5f} {:.5f} \n".format(np.log10(np.sum(gal.star['m'])), np.mean(age), np.mean(gal.star['metal']), np.log10(new_star)))
        
    arr[:3] = theta, phi, alpha
    print(f"new inclination and azimuth: {theta:.3f}, {phi:.3f}")

    if plot_stellar:
        fig, axs = plt.subplots(2,2)
        fig.set_size_inches(6,6)
        fig.suptitle(f"nout={nout}, ID={gid}")
        axs[0,0].hist2d(gal.star["x"], gal.star["y"], bins=256, norm=LogNorm())
        gal.meta.nvec = nvec
        # if it's already along an intrinsic axis, don't bother to rotate it.
        if 1 not in gal.meta.nvec:
            gal.reorient(dest=[0,0,1]) 
        axs[0,1].hist2d(gal.star["x"], gal.star["y"], bins=256, norm=LogNorm())
        axs[1,0].hist2d(gal.star["x"], gal.star["z"], bins=256, norm=LogNorm())
        for ax in axs.ravel():
            ax.set_aspect("equal")
        axs[1,1].text(0.1,0.1,'{:.2f} Msun'.format(np.log10(gcat['m'])), transform=ax.transAxes)
        plt.savefig(dir_out+ f"{gid:05d}_stellar_maps.png", dpi=200)
    
    # update .ski immediately
    fn_template = './template_zubko.ski'
    skifile = sm.SkiFile(fn_template)
    repo = dir_out+ f"{gid:05d}/{nout:05d}/"#faceon_redshift_"
    """
    if inclination + param['incli_off'] < 0:
        off = param['incli_off'] + 180
    elif inclination + param['incli_off'] > 180:
        off = param['incli_off'] - 180
    else:
        off = param['incli_off']
    """
    write_ski(repo, gid, nout, skifile, arr, ski_params)


    print("done") 
