import numpy as np
import load
import utils.sampling as smp
import matplotlib.pyplot as plt


def gm2code(arr, info):
    """
    arr originally is a reference to the object from the outside.
    But later in "arr = arr / info.pboxsize + 0.5", 
    the arr is a new
    """
    return (arr / info.pboxsize + 0.5)# * 0.999783599



def load_gal(nout, idgal):
    s = load.sim.Sim(nout)
    
    # GalaxyMaker dump
    gal = load.rd_GM.rd_gal(nout, idgal, wdir='./')
    
    gal.header['xg'] = gm2code(gal.header['xg'], s.info)
    gal.star['x'] = gm2code(gal.star['x'], s.info)
    gal.star['y'] = gm2code(gal.star['y'], s.info)
    gal.star['z'] = gm2code(gal.star['z'], s.info)
    
    
    # Region.
    radius = 0.5 * max([gal.star['x'].ptp(), gal.star['y'].ptp(), gal.star['z'].ptp()])
    region = smp.set_region(centers=gal.header['xg'], radius=1.5*radius)
    xc, yc, zc = gal.header['xg']
    
    
    # DM particles in the Region.
    s.add_part(ptypes=['dm id pos vel mass'], region=region)
    
    
    # Cell dump
    cell = load.rd_GM.rd_cell(nout, idgal)
    # center on galaxy position
    cell['x'] -= xc 
    cell['y'] -= yc 
    cell['z'] -= zc 
    
    
    kpc_in_cm = 3.0857e21
    msun_in_g = 1.989e33
    gas_mass = cell['var0'] * s.info.unit_d * (cell['dx'] * s.info.boxtokpc * \
               kpc_in_cm)**3 / msun_in_g
    # 'var0' = density. 
    # cell['var'] * s.info.unit_d : density in gram unit (cgs unit system).
    # cell['dx'] * s.info.boxtokpc : cell size in kpc
    # cell['dx'] * s.info.boxtokpc * kpc_in_cm : cell size in cm.
    
    
    star = gal.star
    # center on galaxy position
    star['x'] -= xc 
    star['y'] -= yc 
    star['z'] -= zc 
    star['m'] *=1e11
    
    
    dm = s.part.dm
    ind = np.where(np.square(dm['x'] - xc) + \
                   np.square(dm['y'] - yc) + \
                    np.square(dm['z'] - zc) < np.square(radius))[0]
    dm = dm[ind]
    # center on galaxy position
    dm['x'] -= xc 
    dm['y'] -= yc 
    dm['z'] -= zc 
    dm['m'] *= s.info.msun
    
    
    # Distance
    #rdm  = np.sqrt(dm['x']**2 + dm['y']**2 + dm['z']**2) * s.info.pboxsize * 1e3
    #rst  = np.sqrt(star['x']**2 + star['y']**2 + star['z']**2) * s.info.pboxsize * 1e3
    #rgas = np.sqrt(cell['x']**2 + cell['y']**2 + cell['z']**2) * s.info.pboxsize * 1e3 # in kpc unit.
    
    # Sort by distance
    #rdsort = np.argsort(rdm)
    #rssort = np.argsort(rst)
    #rgsort = np.argsort(rgas)
    
    # cumulative mass sum
    # All mass in Msun unit.
    #cmdm = np.cumsum(dm['m'][rdsort])
    #cmst = np.cumsum(star['m'][rssort])
    #cmgas = np.cumsum(gas_mass[rgsort])
    
    # radial density profile.
    fig, ax = plt.subplots(2)
    #ax.plot(np.log10(rdm[rdsort]), np.log10(cmdm))
    #ax.plot(np.log10(rst[rssort]), np.log10(cmst))
    #ax.plot(np.log10(rgas[rgsort]), np.log10(cmgas))
    #ax.set_xlabel("log(kpc)")
    #ax.set_ylabel("log(Msun)")
    #ax.set_title("Cumulative mass")
    ax[0].hist2d(dm['x'], dm['y'], bins=[100,100])
    ax[1].hist2d(star['x'], star['y'], bins=[100,100])
    plt.show()
    
    
    
