import numpy as np

def get_phi(rsorted,cumsum):
    """
    Binney Ch 2.2
    """
    Grav = 6.67408*1e-11*1.989*1e30/(3.0857*1e16*1e9)

    pi = -Grav*cumsum/rsorted**2*np.append(rsorted[0],rsorted[1:]-rsorted[:-1])
    return np.cumsum(pi[::-1])[::-1]


def get_E(gal, nvec=None, nvec_ys=True, ptype='star', method="Abadi",
                bound_only = True, return_ellip=False, phi_direct=False,
                save_mbp=False, nmbp=100):
    """
    Approximate potential calculation.

    parameters
    ----------
    nvec :
        dd
    nvec_ys :

    ptype : "star"
        Components with which

    save_mbp: false
        Store nmbp most bound star and DM ids.

    nmbp : 100

    """
    # parameters
    pos1, pos2, vel1, vel2 = 'z', 'x', 'vz', 'vx'

    # constants
    G = 6.67384e-11  # m^3 kg^-1 s^-2
    kpc_to_m = 3.08567758e19
    msun_in_kg = 1.9891e30 # kg
    Grav = G * msun_in_kg/(3.0857*1e16*1e9)
    kpc_in_cm = kpc_to_m * 1e2
    msun_in_g = msun_in_kg * 1e3

    info = gal.info

    cell_to_msun = info.unit_d/msun_in_g * kpc_in_cm**3
    cell_mass = gal.cell["rho"]*gal.cell["dx"]**3 * cell_to_msun

    r_s = np.sqrt(np.sum(np.square(gal.star["pos"]), axis=1))
    gal.star = gal.star[np.argsort(r_s)]
    r_s = r_s[np.argsort(r_s)]

    rv = np.cross(gal.star["pos"], gal.star["vel"])

    # Normal vector of the galaxy rotation
    if nvec is None:
        if nvec_ys:
            if gal.star["time"].mean() < 0:
                # it's in conformal time unit.
                from utils.cosmology import Timeconvert
                tc = Timeconvert(info = info)
                gal.star['time'] = tc.time2gyr(gal.star['time'], z_now = info.zred)
                from ..utils.util import replace_field_name
                replace_field_name(gg.star, "time", "age")
            nvec = sum(rv[gal.star['time']<0.01]) # younger than 10Myr
        else:
            R90m = get_radius(gal.star,gal.star['m'],0.9)
            nvec = sum(rv[r_s<R90m])

    # Calculating boundness requires total mass inside a radius.
    # -> DM, Cell are also needed.
    if hasattr(gal, "cell"):
        m_g = cell_mass
        m_all = np.concatenate((gal.star['m'],
                                m_g,
                                gal.dm['m']))
        r_all = np.concatenate((r_s,
                                np.sqrt(np.sum(np.square(gal.cell["pos"]), axis=1)),
                                np.sqrt(np.sum(np.square(gal.dm  ["pos"]), axis=1))))
    else:
        m_all = np.concatenate((gal.star['m'],
                                gal.dm['m']))
        r_all = np.concatenate((np.sqrt(np.sum(np.square(gal.star["pos"]), axis=1)),
                                np.sqrt(np.sum(np.square(gal.dm  ["pos"]), axis=1))))

    i_sorted = np.argsort(r_all)
    r_all_sorted = r_all[i_sorted]
    m_enc = np.cumsum(m_all[i_sorted])

    # First nstar indices are stars.

    # Exclude the star at the center.
    x  = gal.star[pos1]
    y  = gal.star[pos2]
    vx = gal.star[vel1]
    vy = gal.star[vel2]
    m  = gal.star['m']

    #boxtokpc = gal.info.pboxsize * 1000
    i_star_dist = np.searchsorted(r_all_sorted, r_s)

    if method=="scannapieco":
        v_circ = np.sqrt(G * msun_in_kg * m_enc[i_star_dist]/
                         (kpc_to_m * r_s)) * 1e-3 # m/s to in km/s
        j_circ = r_s * v_circ
        # Finally, r in code unit, v in km/s
        j_phi = np.inner(rv, [0,1,0])
        # Equivalent, but below is probably faster.
        #j_phi = (x * vy - y * vx) # * boxtokpc omitted.
        gal.E = j_phi / j_circ
    elif method=="Abadi":
        #r_s_sorted = r_s[1:]

        RV = np.inner(rv,nvec)#/np.sqrt(sum(J*J))
        cos_alp = RV/np.sqrt((rv*rv).sum(axis=1))

        st_ind = np.searchsorted(r_all_sorted, r_s)
        M_incl_star = m_enc[st_ind][1:]

        jz = RV[1:]
        vcir = np.sqrt(Grav*M_incl_star/r_s[1:])
        jcir = r_s[1:]*vcir
        e = jz/jcir
        cos = cos_alp[1:]

        if phi_direct:
            from analysis import pot_diret
            Grav = 6.67408*1e-11*1.989*1e30/(3.0857*1e16*1e9)
            #md,xd,yd,zd,mg,xg,yg,zg,ms,xs,ys,zs,eps,nd,ng,ns
            phi = pot_diret.star_potential(gal.dm['m'],
                                     gal.dm['x'],
                                     gal.dm['y'],
                                     gal.dm['z'],
                                     cell_mass,
                                     gal.cell["x"],
                                     gal.cell["y"],
                                     gal.cell["z"],
                                     gal.star["m"],
                                     gal.star["x"],
                                     gal.star["y"],
                                     gal.star["z"],
                                     np.ones_like(gal.star["m"])*gal.cell["dx"].min()*0.3,
                                     len(gal.dm), len(gal.cell), len(gal.star))
            phi *= -Grav
        else:
            phi = get_phi(r_s[1:], M_incl_star)

        specificE = phi + 0.5*(gal.star["vel"][:,0]**2
                               +gal.star["vel"][:,1]**2
                               +gal.star["vel"][:,2]**2)[1:]

        if save_mbp:
            gal.meta.mbp_star = gal.star["id"][np.argsort(specificE)[:nmbp]]

            r_dm = np.sqrt(np.sum(np.square(gal.dm["pos"]), axis=1))
            i_dm_sort = np.argsort(r_dm)
            gal.dm = gal.dm[i_dm_sort]
            r_dm = r_dm[i_dm_sort]

            dm_ind = np.searchsorted(r_all_sorted, r_dm)
            M_incl_dm = m_enc[dm_ind]
            phi_DM = get_phi(r_dm, M_incl_dm)
            specificE_DM = phi_DM + 0.5*(gal.dm["vel"][:,0]**2
                                        +gal.dm["vel"][:,1]**2
                                        +gal.dm["vel"][:,2]**2)
            gal.meta.mbp_DM = gal.dm["id"][np.argsort(specificE_DM)[:nmbp]]

        Ecir = 0.5*Grav*M_incl_star/r_s[1:] + phi

        i_sort_Ecir = np.argsort(Ecir)
        jE = jcir[i_sort_Ecir][np.digitize(specificE,Ecir[i_sort_Ecir])-1]

        try:
            gal.star["ellip"][1:] = jz/jE
            gal.star["ellip"][0] = 0
        except:
            print("[get_E] ellip field is not available")
            gal.E = jz/jE

        gal.e = e
        gal.cos = cos
        gal.vcir = vcir
        gal.id_closest_star = gal.star["id"][0]
