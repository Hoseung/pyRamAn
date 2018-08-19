import numpy as np
from utils.sampling import Region

def mk_gal(gal,
            save=False,
            verbose=False,
            mstar_min=1e9,
            den_lim=1e6,
            den_lim2=5e6,
            rmin = -1,
            Rgal_to_reff=5.0,
            method_com="catalog",
            method_cov="catalog",
            method_member="Reff",
            follow_bp=False,
            unit_conversion="code",
            convert_time=False):
    """
        Determine if this is a legitimate galxy. re-center components.

        This routine consists of three parts:
        1) decide if the system is dense enough.
        2) determine member components
        3) re-center components and do simple calculations (com, cov, total mass, ...)

        Raw star/DM/gas data are only rough representations of a galaxy.
        But not much have to be thrown away, either.

        Parameters
        ----------

        Rgal_to_reff:
            Galaxy radius = Reff * Rgal_to_reff.
            By default, Rgal = 5*Reff.
            (E's have smaller R_t_r, L's have larger.)
        save: False,
        verbose:False,
        mstar_min:1e9,
        den_lim:1e6,
        den_lim2:5e6,
        rmin : -1,
        Rgal_to_reff:5.0,
        method_com:"catalog",
        method_cov:"catalog",
        method_member:"Reff",
        follow_bp:False,
        unit_conversion:"code"
        convert_time:
            default = False
            because it is more efficient to create an timecoverter instance once (load table) and
            convert multiple galaxies altogether.

        Notes
        -----
        1. Since now (as of 2016.03) galaxy calculation is based on the GalaxyMaker results,
        let's just believe and minimize redundant processes such as determining the center of mass,
        system velocity, and checking the presence of (unwanted) substructure.

        2. Assuming all catalog in the code units.

        3. dr, bin size in determining the radial profile cut scales with exponential factor.
        Without scaling, 0.5kpc bin size is too large for ~1kpc galaxies at high-z.

    """
    # Need halo information
    assert (gal.gcat is not None), ("Need a catalog,"
    "use Galaxy.set_halo() and provide x,y,z,vx,vy,vz,r,rvir, at least"
    "Units are.. ?? ")

    pbx = gal.info.pboxsize

    # galaxy center from GalaxyMaker. - good enough.
    xc = gal.gcat["x"]
    yc = gal.gcat["y"]
    zc = gal.gcat["z"]
    if verbose:
        print("Galaxy center : {} {} {} using {}".format(xc, yc, zc, method_com))

    vxc = gal.gcat["vx"]
    vyc = gal.gcat["vy"]
    vzc = gal.gcat["vz"]
    if verbose:
        print("Velocity center : {} {} {} using {}".format(vxc,vyc,vzc,method_cov))

    star = gal.star
    # re-center position first.
    #if star is not None:
    #    star["x"] = (star["x"] - xc)*1e3
    #    star["y"] = (star["y"] - yc)*1e3
    #    star["z"] = (star["z"] - zc)*1e3
    #    star["m"] *= 1e11#gal.info.msun
    if verbose: print("star x", gal.star["x"])

    dm = gal.dm
    if dm is not None:
        gal._has_dm = True
        dm["x"] = (dm["x"] - xc)*1e3
        dm["y"] = (dm["y"] - yc)*1e3
        dm["z"] = (dm["z"] - zc)*1e3
        dm["m"] *= gal.info.msun

    cell = gal.cell
    # Don't convert cell units here,
    # COPY only relevant cells and then modify them.
    if cell is not None:
        gal.__has_cell = True

    assert (gal._has_star or gal._has_dm or gal._has_cell), ("At least"
    "one of three(star, dm, gas) component is needed")

    # all tests passed.
    if verbose:
        print("Making a galaxy:", gal.meta.id)
        print("SAVE:", save)
        print("Halo size:", gal.gcat['rvir'])


    rgal_tmp = min([gal.gcat['r'] * 1e3, 30]) # gcat["rvir"] in kpc
    if verbose:
        print("Rgal_tmp", rgal_tmp)
        print("gal.debug",gal.debug)
    dense_enough = radial_profile_cut(gal, star['x'], star['y'], star['m'],
                         den_lim=den_lim, den_lim2=den_lim2,
                         mag_lim=25,
                         nbins=int(rgal_tmp/0.5),
                         dr=0.5 * gal.info.aexp,
                         rmax=rgal_tmp,
                         debug=gal.debug)
    if not dense_enough:
        print("Not dense enough")
        return False

    if method_com=="catalog":
        gal.meta.xc, gal.meta.yc, gal.meta.zc = gal.header["xg"]

    dd = (np.square(star['x']) +
          np.square(star['y']) +
          np.square(star['z']))

    if method_cov=="close_member":
        i_close = np.argsort(dd)[:int(len(star))] # half close members
        gal.meta.vxc = np.average(star["vx"][i_close])
        gal.meta.vyc = np.average(star["vy"][i_close])
        gal.meta.vzc = np.average(star["vz"][i_close])
        print("method_COV = close_member")
    elif method_cov=="catalog":
        gal.meta.vxc, gal.meta.vyc, gal.meta.vzc = gal.header["vg"]

    # Membership
    ind = np.where(dd < gal.meta.rgal**2)[0]# in kpc unit

    gal.star = star[ind]
    gal.star["vx"] -= gal.meta.vxc
    gal.star["vy"] -= gal.meta.vyc
    gal.star["vz"] -= gal.meta.vzc

    if gal.debug:
        print('[galaxy.Galaxy.mk_gal] mimax vx :',
              min(gal.star['vx']),
              max(gal.star['vx']))

    gal.meta.nstar = len(ind)
    gal.meta.mstar = sum(gal.star['m'])

    if gal.meta.mstar < mstar_min:
        print("Not enough stars: {:.2e} Msun".format(gal.meta.mstar))
        print("{} Aborting... \n".format(len(gal.star['m'])))
        gal.meta.star = False
        return False

    gal.meta.Rgal_to_reff = gal.meta.rgal / gal.meta.reff
    # should not surpass rr_tmp, over where another galaxy might be.

    # Test if the outer annulus has significant amount of stars
    # -> it shouldn't.
    if star is not None:
        nstar_tot = len(star['x'][ind])
        if verbose: print("nstar tot:", nstar_tot)
        if verbose: print("Store stellar particle")

        if 'time' in gal.star.dtype.names and convert_time:
            from utils.cosmology import Timeconvert
            tc = Timeconvert(gal.info)
            gal.star['time'] = tc.time2gyr(gal.star['time'],
                                 z_now = gal.info.zred)

    # VERY arbitrary..
    rgal_tmp = gal.meta.Rgal_to_reff *gal.meta.reff

    #print(".........", gal.star['m'][100:120], gal.mstar)

    gal.region = Region(xc=gal.meta.xc,
                        yc=gal.meta.yc,
                        zc=gal.meta.zc,
                        radius = gal.meta.rgal)

    if gal.debug:
        print('[galaxy.Galaxy.mk_gal] meta.v[x,y,z]c',
              gal.meta.vxc, gal.meta.vyc, gal.meta.vzc)
        print('[galaxy.Galaxy.mk_gal] mima vx 2',
              min(gal.star['vx']), max(gal.star['vx']))

    if dm is not None:
        if method_member == "Reff":
            idm = np.where( np.square(dm["x"] - gal.meta.xc) +
                            np.square(dm["y"] - gal.meta.yc) +
                            np.square(dm["z"] - gal.meta.zc) <= np.square(rgal_tmp))[0]
        elif method_member == "v200":
        # Although the velocity is redefined later,
        # particle membership is fixed at this point.
            idm = np.where( np.square(dm["vx"] - gal.meta.vxc / gal.info.kms)+
                            np.square(dm["vy"] - gal.meta.vyc / gal.info.kms)+
                            np.square(dm["vz"] - gal.meta.vzc / gal.info.kms) <= np.square(200**2))[0]
        gal._convert_unit("dm", unit_conversion)

    if cell is not None:
        dtype_cell = [('x', '<f8'), ('y', '<f8'), ('z', '<f8'),
                      ('dx', '<f8'), ('rho', '<f8'),
                      ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),
                      ('temp', '<f8')]#, ('metal', '<f8')]
        if "var5" in cell.dtype.names:
            if len(cell.dtype.names) < 12:
                dtype_cell.append(("metal", "<f8"))
            else:
                print("[mk_gal] Warning...")
                print("[mk_gal] Don't know what to do with all the hydro-variables:")
                print("[mk_gal] ",cell.dtype)
                print("[mk_gal] Ignoring anyting after the temperature field.")

        if "cpu" in cell.dtype.names:
            dtype_cell.append(('cpu', '<f8'))

        if verbose: print("Cell is NOT none")
        icell = np.where(np.square(cell["x"] - (xc/pbx + 0.5)) +
                         np.square(cell["y"] - (yc/pbx + 0.5)) +
                         np.square(cell["z"] - (zc/pbx + 0.5)) <= np.square(rgal_tmp))[0]

        #gal._add_cell(cell, icell)
        #gal._convert_unit("cell", unit_conversion)

        gal.cell = np.recarray(len(icell), dtype=dtype_cell)
        gal.cell['x'] = (cell['x'][icell] - 0.5) * pbx * 1e3 - xc
        gal.cell['y'] = (cell['y'][icell] - 0.5) * pbx * 1e3 - xc
        gal.cell['z'] = (cell['z'][icell] - 0.5) * pbx * 1e3 - xc
        gal.cell['dx'] = cell['dx'][icell] * pbx * 1000
        gal.cell['rho'] = cell['var0'][icell]
        gal.cell['vx'] = cell['var1'][icell] * gal.info.kms - gal.meta.vxc
        gal.cell['vy'] = cell['var2'][icell] * gal.info.kms - gal.meta.vyc
        gal.cell['vz'] = cell['var3'][icell] * gal.info.kms - gal.meta.vzc
        gal.cell['temp'] = cell['var4'][icell]
        if "var5" in cell.dtype.names:
            gal.cell['metal'] = cell['var5'][icell]
        if "cpu" in cell.dtype.names:
            gal.cell['cpu'] = cell['cpu'][icell]

        print("cell x, final", cell["x"])
        gal.cal_mgas()

    # Some more sophistications.
    """
    print("Rgal = 4 * Reff = ", rgal_tmp * gal.info.pboxsize * 1000)

        # Save sink particle as a BH, not cloud particles.

    """
    return True


def extract_cold_gas(gg, rmax = 180, dr = 5):
    """
        Measure radial profile and returns indices of cells inside r_min,
        where r_min is the local minima of radial MASS profile.
        -> should I use density profile instead?
    """
    from scipy.signal import argrelmin
    # radial profile.
    if not hasattr(gg,"cell"):
        print("No cell found")
        return

    cold_cell = gg.cell[rho_t_cut(gg.cell, gg.info)]


    rr = np.sqrt(np.square(cold_cell["x"])+\
                 np.square(cold_cell["y"])+\
                 np.square(cold_cell["z"]))

    i_sort = np.argsort(rr)
    r_sorted = rr[i_sort]
    mm = cold_cell["dx"]**3 * cold_cell["var0"]
    m_sorted = mm[i_sort]

    rmax = min([np.max(rr), rmax])

    # Note 1.
    # Depends on the cell resolution. How about 8 * dx_min?
    # Larger dx will count in small satellites,
    # while smaller dx will make the measurement sensitive to density fluctuations.
    nbins= int(rmax/dr)

    frequency, bins = np.histogram(r_sorted, bins = nbins, range=[0, rmax])
    bin_centers = bins[:-1] + 0.5 * dr # remove the rightmost boundary.

    m_radial = np.zeros(nbins)
    ibins = np.concatenate((np.zeros(1,dtype=int), np.cumsum(frequency)))

    for i in range(nbins):
        m_radial[i] = np.sum(m_sorted[ibins[i]:ibins[i+1]])
        # Check stellar surface density
        sig_at_r = m_radial[i]/(2 * np.pi * bin_centers[i] * dr)

    # Find local minimum
    # If there is flat zeros, take the first zero.
    # If not, use scipy.argrelmin
    i_zero = np.argmax(m_radial==0)
    if i_zero > 0:
        ind_min = i_zero -1
    else:
        ind_min= argrelmin(m_radial)[0] -1 # 1D array for 1D input.
        ind_min = ind_min[np.argmax(ind_min * dr > rmin)]* dr

    # Note 2.
    # If the minimum is farther than rmin=10kpc,
    # I assume that is correct.
    gg.cell = cold_cell[rr < bin_centers[ind_min]]
    gg.mgas_cold = np.sum(gg.cell["var0"]*gg.cell["dx"]**3)
    gg.cold_gas_profile = dict(profile=m_radial[:ind_min],dr=dr)

def rho_t_cut(cell, info, clip_sigma=0):
    """
        Extract galactic cold gas following Torrey+12 criterion.
        Assume cells in the original (code) unit.
    """
    # Var0 in Msun h^2 kpc^-3 unit.
    kpc_in_cm = 3.08567758e21
    msun_in_g = 1.99e33
    gcc2this_unit = kpc_in_cm**3/msun_in_g
    if clip_sigma > 0:
        pass
        #Do sigma clipping..
    return np.log10(cell["var4"]/cell["var0"]*info.unit_T2) < 6 + 0.25*np.log10((cell["var0"]*info.unit_d)*gcc2this_unit*1e-10) #


def radial_profile_cut(gal, xx, yy, mm,
                       den_lim=1e6, den_lim2=5e6,
                       mag_lim=25, nbins=100, rmax=20, dr=0.5,
                       debug=False):
    """
    System velocity determined as np.average(vx[i_close]) sometimes fail,
    which may indicate that this function fails to extract reliable member stars.
    This occurs more frequently with high-z or high resolution data.
    Todo
    ----
        Adaptive member determination over varying resolution and redshift.


    """

    rr = np.sqrt(np.square(xx) + np.square(yy))# in kpc unit
    if debug:
        print("min(rr) {}\n max(rr){}\n min(xx){}\n max(xx){}".format(
                               min(rr), max(rr), min(xx), max(xx)))

    # Mass weight.
    i_sort = np.argsort(rr)
    r_sorted = rr[i_sort]
    m_sorted = mm[i_sort]

    rmax = min([np.max(rr), rmax])
    nbins = int(rmax/dr)

    if nbins < 3:
        print("Too small size \n # of stars:", len(rr))
        return False

    frequency, bins = np.histogram(r_sorted, bins = nbins, range=[0, rmax])
    bin_centers = bins[:-1] + 0.5 * dr # remove the rightmost boundary.

    m_radial = np.zeros(nbins)
    ibins = np.concatenate((np.zeros(1,dtype=int), np.cumsum(frequency)))

    i_r_cut1 = nbins -1 # Maximum value
    # on rare occasions, a galaxy's stellar surface density
    # never crosses the density limit. Then i_r_cut1 = last index.
    for i in range(nbins):
        m_radial[i] = np.sum(m_sorted[ibins[i]:ibins[i+1]])
        # Check stellar surface density
        sig_at_r = m_radial[i]/(2 * np.pi * bin_centers[i] * dr)
        if debug:
            print(sig_at_r, den_lim)
        if sig_at_r < den_lim:
            i_r_cut1 = i-1
            break
    #i_r_cut2= np.argmax(m_radial/(2 * np.pi * bin_centers * dr) < den_lim2)
    # If for some reason central region is less dense,
    # profile can end at the first index.
    # Instead coming from backward, search for the point the opposite condition satisfied.
    if debug:
        print(rmax, nbins)
        print("frequency", frequency)
        print("bins", bins)
        print("ibins", ibins)
        print("bin centers", bin_centers)
        print("m_radial", m_radial)

    den_radial_inverse = m_radial[::-1]/(2 * np.pi * bin_centers[::-1] * dr)
    if debug: print("den_radial_inverse", den_radial_inverse)
    if max(den_radial_inverse) < 2 * den_lim2:
        np.set_printoptions(precision=3)
        print("radial density :",den_radial_inverse)
        return False
    i_r_cut2=len(m_radial) - np.argmax(den_radial_inverse > den_lim2) -1
    if debug:
        print("[galaxy.Galaxy.radial_profile_cut] m_radial \n", m_radial)
        print("[galaxy.Galaxy.radial_profile_cut] den_radial_inverse \n", den_radial_inverse)
        print("[galaxy.Galaxy.radial_profile_cut] i_r_cut2", i_r_cut2)

    mtot2 = sum(m_radial[:i_r_cut2])
    mtot1 = sum(m_radial[:i_r_cut1])
    i_reff2 = np.argmax(np.cumsum(m_sorted) > (0.5*mtot2))
    i_reff1 = np.argmax(np.cumsum(m_sorted) > (0.5*mtot1))
    gal.meta.reff2 = r_sorted[i_reff2]
    gal.meta.reff  = r_sorted[i_reff1]
    gal.meta.rgal2 = max([bin_centers[i_r_cut2],4*gal.meta.reff2])
    gal.meta.rgal  = max([bin_centers[i_r_cut1],4*gal.meta.reff])#bin_centers[i_r_cut1]

    #       It is not wrong for BCGs to have very large Reff(~50kpc).
    #       But referring the average velocity of stellar particles inside 50kpc
    #       as the system velocity is WRONG.
    #       If 1Reff is huge, try smaller aperture when measuring the system velocity.

    if debug:
        print("[galaxy.Galaxy.radial_profile_cut] mtot, mtot2", mtot1, mtot2)

    i_close = i_sort[:np.argmax(np.cumsum(m_sorted) > (0.2*mtot2))] # 20% closest particles

    return True


def get_normal_vec(stars):
    #vec_rot = np.cross(stars["pos"],stars["vel"])
    vec_rot = np.cross(np.vstack((stars["x"],stars["y"],stars["z"])),
                       np.vstack((stars["vx"],stars["vy"],stars["vz"])))
    Ln = vec_rot.sum(axis=0)
    return Ln / np.linalg.norm(Ln)
