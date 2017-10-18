import load
import numpy as np
import pickle
import matplotlib.pyplot as plt

read=True


class Meta():
    pass

def radial_profile_cut(meta, xx, yy, mm, vx, vy, vz,
                       den_lim=1e6, den_lim2=5e6,
                       mag_lim=25, nbins=100, rmax=50, dr=0.5):
    # 2D photometry. (if rotated towards +y, then use x and z)
    # now assuming +z alignment. 
    rr = np.sqrt(np.square(xx) + np.square(yy))# in kpc unit

    # Account for weights.
    i_sort = np.argsort(rr)
    r_sorted = rr[i_sort]
    m_sorted = mm[i_sort]
 
    rmax = np.max(rr)
    nbins = int(rmax/dr)

    frequency, bins = np.histogram(r_sorted, bins = nbins, range=[0,rmax])
    bin_centers = bins[:-1] + 0.5 * dr # remove the rightmost boundary.
 
    m_radial = np.zeros(nbins)
    ibins = np.concatenate((np.zeros(1), np.cumsum(frequency)))
    for i in range(nbins):
        m_radial[i] = np.sum(m_sorted[ibins[i]:ibins[i+1]])
        if (m_radial[i]/(2 * np.pi * bin_centers[i] * dr)) < den_lim:
            i_r_cut1 = i-1
            break
 
    i_r_cut2= np.argmax(m_radial/(2 * np.pi * bin_centers * dr) < den_lim2)

    mtot2 = sum(m_radial[:i_r_cut2])
    mtot1 = sum(m_radial[:i_r_cut1])
    i_reff2 = np.argmax(np.cumsum(m_sorted) > (0.5*mtot2))
    i_reff1 = np.argmax(np.cumsum(m_sorted) > (0.5*mtot1))
    meta.reff2 = r_sorted[i_reff2]
    meta.reff  = r_sorted[i_reff1]
    meta.rgal2 = max([bin_centers[i_r_cut2],4*meta.reff2])
    meta.rgal  = max([bin_centers[i_r_cut1],4*meta.reff])#bin_centers[i_r_cut1]

    # velocity center
    # It is not wrong for BCGs to have very large Reff(~50kpc). 
    # But referring the average velocity of stellar particles inside 50kpc 
    # as the system velocity is WRONG.
    # If 1Reff is huge, try smaller aperture when measuring the system velocity.

    i_close = i_sort[:np.argmax(np.cumsum(m_sorted) > (0.1*mtot2))] # 10% closest particles
    meta.vxc = np.average(vx[i_close])
    meta.vyc = np.average(vy[i_close])
    meta.vzc = np.average(vz[i_close])


def rp(gal, meta):
    radial_profile_cut(meta, gal.data['pos'][:,0]
                           , gal.data['pos'][:,1]
                           , gal.data['pos'][:,2]
                           , gal.data['vel'][:,0]
                           , gal.data['vel'][:,1]
                           , gal.data['vel'][:,2])


nout = 187
wdir = './29176/'

if read:
    sig3d=[]
    mass=[]
    reff=[]
    import tree.halomodule as hmo
    gg = hmo.Halo(base=wdir, nout=nout, is_gal=True, load=True)
    ngal = len(gg.data)
    print("Total {} galaxies in {}".format(ngal, wdir))
    for igal in range(1,ngal + 1):
        gal = load.rd_GM.rd_gal(nout, igal, wdir=wdir)
        sig3d.append((np.std(gal.data['vel'][:,0])+
                      np.std(gal.data['vel'][:,0])+
                      np.std(gal.data['vel'][:,0]))/np.sqrt(3))
        mass.append(gal.header['mgal'])
 
        meta = Meta()
        rp(gal, meta)
        reff.append(meta.reff)

    mass = np.array(mass)
    sig3d = np.array(sig3d)
    reff = np.array(reff)
    pickle.dump((mass, sig3d, reff),open("mass_sig_r.pickle", 'wb'))
else:
    mass, sig3d, reff = pickle.load(open("mass_sig_r.pickle", 'rb'))

fig, ax = plt.subplots()
G = 6.67408e-11 #[m3 kg-1 s-2]

sig = sig3d * 1e3 # [m/s]
msun = 1.989e30 #[kg]
mass = mass * msun
kpc_to_m = 3.0857e16 #[m]
reff = reff * kpc_to_m

ax.scatter(G * mass, np.array(sig)**2 * np.array(reff))
ax.set_xscale('log')
ax.set_yscale('log')

#plt.show()
plt.savefig("Mass_vs_sigR.pdf")


# Roughly, G*M = R*sig^2. (virialized)


