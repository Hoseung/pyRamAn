import numpy as np

def get_phi(rsorted,cumsum):
    """
    Binney Ch 2.2
    """
    Grav = 6.67408*1e-11*1.989*1e30/(3.0857*1e16*1e9)

    pi = -Grav*cumsum/rsorted**2*np.append(rsorted[0],rsorted[1:]-rsorted[:-1])
    return np.cumsum(pi[::-1])[::-1]


def get_E(self, nvec=None, nvec_ys=True, ptype='star', method="Abadi",
                bound_only = True, return_ellip=False):

    # parameters
    pos1, pos2, vel1, vel2 = 'z', 'x', 'vz', 'vx'

    # constants
    G = 6.67384e-11  # m^3 kg^-1 s^-2
    kpc_to_m = 3.08567758e19
    msun_in_kg = 1.9891e30 # kg
    Grav = G * msun_in_kg/(3.0857*1e16*1e9)
    kpc_in_cm = kpc_to_m * 1e2
    msun_in_g = msun_in_kg * 1e3

    info = self.info

    cell_to_msun = info.unit_d/msun_in_g * kpc_in_cm**3
    cell_mass = self.cell["rho"]*self.cell["dx"]**3 * cell_to_msun

    r_s = np.sqrt(np.sum(np.square(self.star["pos"]), axis=1))
    if False:
        i_ok = np.zeros(len(r_s), dtype="bool")
        i_ok[:] = True
        i_ok[np.argmin(r_s)] = False

        star = self.star[i_ok]
        r_s = r_s[i_ok]
    else:
        star = self.star

    rv = np.cross(star["pos"], star["vel"])
    if nvec is None:
        if nvec_ys:
            if star["time"].mean() < 0:
                # it's in conformal time unit.
                from utils.cosmology import Timeconvert
                tc = Timeconvert(info = info)
                star['time'] = tc.time2gyr(star['time'], z_now = info.zred)
            nvec = sum(rv[star['time']<0.01]) # younger than 10Myr
        else:
            R90m = get_radius(star,star['m'],0.9)
            nvec = sum(rv[r_s<R90m])


    # Calculating boundness requires total mass inside a radius.
    # -> DM, Cell are also needed.
    if hasattr(self, "cell"):
        m_g = cell_mass
        m_all = np.concatenate((self.star['m'],
                                m_g,
                                self.dm['m']))
        r_all = np.concatenate((r_s,
                                np.sqrt(np.sum(np.square(self.cell["pos"]), axis=1)),
                                np.sqrt(np.sum(np.square(self.dm  ["pos"]), axis=1))))
    else:
        m_all = np.concatenate((self.star['m'],
                                self.dm['m']))
        r_all = np.concatenate((np.sqrt(np.sum(np.square(self.star["pos"]), axis=1)),
                                np.sqrt(np.sum(np.square(self.dm  ["pos"]), axis=1))))

    i_sorted = np.argsort(r_all)
    r_all = r_all[i_sorted]
    m_enc = np.cumsum(m_all[i_sorted])

    # First nstar indices are stars.

    # Exclude the star at the center.
    bound_only = False

    if bound_only:
        i_star = i_sorted[self.bound_ptcl]
        x  = self.star[pos1][self.bound_ptcl]
        y  = self.star[pos2][self.bound_ptcl]
        vx = self.star[vel1][self.bound_ptcl]
        vy = self.star[vel2][self.bound_ptcl]
        m  = self.star['m'] [self.bound_ptcl]
    else:
        #i_star = i_sorted[:len(r_s)]
        x  = star[pos1]
        y  = star[pos2]
        vx = star[vel1]
        vy = star[vel2]
        m  = star['m']


    #boxtokpc = self.info.pboxsize * 1000
    i_star_dist = np.searchsorted(r_all, r_s)

    if method=="scannapieco":
        v_circ = np.sqrt(G * msun_in_kg * m_enc[i_star_dist]/
                         (kpc_to_m * r_s)) * 1e-3 # m/s to in km/s
        j_circ = r_s * v_circ
        # Finally, r in code unit, v in km/s
        j_phi = np.inner(rv, [0,1,0])
        # Equivalent, but below is probably faster.
        #j_phi = (x * vy - y * vx) # * boxtokpc omitted.
        self.E = j_phi / j_circ
    elif method=="Abadi":
        rssort = np.argsort(r_s)
        r_s_sorted = r_s[rssort][1:]

        RV = np.inner(rv,nvec)#/np.sqrt(sum(J*J))
        cos_alp = RV/np.sqrt((rv*rv).sum(axis=1))

        st_ind = np.searchsorted(r_all, r_s[rssort])

        jz = RV[rssort][1:]
        M = m_enc[st_ind][1:]

        vcir = np.sqrt(Grav*M/r_s[rssort][1:])
        jcir = r_s[rssort][1:]*vcir
        e = jz/jcir
        cos = cos_alp[rssort][1:]

        phi = get_phi(r_s_sorted, M)
        specificE = phi + 0.5*(star["vel"][:,0]**2+star["vel"][:,1]**2
                               +star["vel"][:,2]**2)[rssort][1:]
        Ecir = 0.5*Grav*M/r_s_sorted + phi

        i_sort_Ecir = np.argsort(Ecir)
        jE = jcir[i_sort_Ecir][np.digitize(specificE,Ecir[i_sort_Ecir])-1]

        try:
            self.star["ellip"][rssort[1:]] = jz/jE
            self.star["ellip"][rssort[0]] = 0
        except:
            print("[get_E] ellip field is not available")
            self.E = jz/jE

        #gal.sorted = star[rssort][1:]
        self.e = e
        self.cos = cos
        self.vcir = vcir
        self.ind_closest_star = rssort[0]
