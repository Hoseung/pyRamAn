import numpy as np

# 1. Gas properties
def get_gas_all(gg, info, dr=5, rmax=200, density_ratio=1e-3):
    #get_cold_cell(gg, info, dr, rmax, density_ratio)
    #mgas_tot, mgas_cold, Ln_gas = get_gas_properties(gg, info)
    gg.meta.gas_results = get_gas_properties(gg, info)


def get_axis_ratio(gg):
    cov_mat = np.cov((gg.star["x"], gg.star["y"], gg.star["z"]))
    eig_val_cov, eig_vec_cov = np.linalg.eig(cov_mat)
    gg.meta.abc_eig_vec = np.sort(eig_val_cov)


def get_cold_cell(gg, info, dr=5, rmax=200, density_ratio=1e-3, check_bound=False):
    """
        returns True if the galaxy has cold gas.
    """
    cold_cell = gg.cell[rho_t_cut(gg.cell, info)]
    if len(cold_cell) < 1:
        # this galaxy has no cold gas. That's possible.
        return False
    else:
        #print("Numbe of cold cells", len(cold_cell))
        if check_bound:
            i_dense = ind_dense(cold_cell, dr=dr, rmax=rmax)
            if len(i_dense)/len(cold_cell) < density_ratio:
                print("Warning.. only a tiny fraction of cold gas is bound to the galaxy?")
            gg.cell = cold_cell[i_dense]
        else:
            gg.cell = cold_cell
        return True

def rho_t_cut(cell, info, lose_cut=False):
    """
        Extract galactic cold gas following Torrey+12 criterion.
        Assume cells in the original (code) unit.
    """
    # Var0 in Msun h^2 kpc^-3 unit.
    kpc_in_cm = 3.08567758e21
    msun_in_g = 1.99e33
    gcc2this_unit = kpc_in_cm**3/msun_in_g

    return np.log10(cell["var4"]/cell["var0"]*info.unit_T2) < 6 + 0.25*np.log10((cell["var0"]*info.unit_d)*gcc2this_unit*1e-10)#


def ind_dense(cell, rmax = 200, dr = 5, rmin=1):
    """
        Measure radial profile and returns indices of cells inside r_min,
        where r_min is the local minima of radial MASS profile.
        -> should I use density profile instead?
    """
    from scipy.signal import argrelmin
    # radial profile.
    rr = np.sqrt(np.square(cell["x"])+\
                 np.square(cell["y"])+\
                 np.square(cell["z"]))

    i_sort = np.argsort(rr)
    r_sorted = rr[i_sort]
    mm = cell["dx"]**3 * cell["var0"]
    m_sorted = mm[i_sort]
    rmax = max([10, min([np.max(rr), rmax])])
    #print(rmax)

    #print("rmax now", rmax)
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
    # 1. If there is flat zeros, take the first zero.
    # If not, use scipy.argrelmin
    i_zero = np.argmax(m_radial==0)
    if i_zero > 0:
        ind_min = i_zero -1
    else:
        try:
            ind_min= min(argrelmin(m_radial)[0]) -1 # 1D array for 1D input.
        except:
            ind_min = -1
        #ind_min = ind_min[np.argmax(ind_min * dr > rmin)]* dr

    # Note 2.
    # If the minimum is farther than rmin=10kpc,
    # I assume that is correct.
    return rr < bin_centers[ind_min]


def get_gas_properties(gg, info):
    """
        Calculates cold gas mass, angular

        Length must be in kpc unit, velocities in km/s unit.
        But densities... not sure yet.


    """
    vrho = "var0"
    vdx = "dx"
    vvx = "var1"
    vvy = "var2"
    vvz = "var3"

    # total gas mass before extracting cold gas.
    gg.meta.gas_results["mgas_tot"] = np.sum(gg.cell[vrho]*gg.cell[vdx]**3)

    # Leave only cold gas and measure properties.
    has_cold_gas = get_cold_cell(gg, info)
    if has_cold_gas:
        gg.meta.gas_results["mgas_cold"] = np.sum(gg.cell[vrho]*gg.cell[vdx]**3)
        vec_rot = np.cross(np.stack((gg.cell["x"],gg.cell["y"],gg.cell["z"])).T,
                           np.stack((gg.cell[vvx],
                                     gg.cell[vvy],
                                     gg.cell[vvz])).T)

        gg.meta.gas_results["Ln_gas"]=(vec_rot.T * (gg.cell[vdx]**3 *gg.cell[vrho])).sum(axis=1)
    else:
        gg.meta.gas_results["mgas_cold"]= 0
        gg.meta.gas_results["Ln_gas"] = (-1,-1,-1)


# 2. SFR
def get_sfr_all(gg, sfr_dts=0.1,
                hist_dt=0.1, hist_tmin=0, hist_tmax=None):
    h = get_age_hist(gg, dt=hist_dt, tmin=hist_tmin, tmax=hist_tmax)
    sfrs = get_sfr(gg, sfr_dts, age_hist=h, dt_age_hist=hist_dt)
    area = get_galaxy_area(gg)

    gg.meta.sfr_results["sfr_dts"]=sfr_dts
    gg.meta.sfr_results["sfrs"]=sfrs/area
    gg.meta.sfr_results["area"]=area
    gg.meta.sfr_results["hist"]=h/area
    gg.meta.sfr_results["hist_dt"]=hist_dt
    gg.meta.sfr_results["hist_tmin"]=hist_tmin
    gg.meta.sfr_results["hist_tmax"]=hist_tmax


def get_age_hist(gg,
                dt = 0.1,
                tmin=0,
                tmax=None):
    """
        Returns results in the result_sfr format.
        I.e., dt, tmin, tmax, and the histogram in Msun * yr-1  unit.

        Note
        ----
        the result container namedtuple is NOT generated inside the funtion.
        the result container subclass is defined in the outer function once
        and each instance is generated for each galaxy.
        This suppresses defining the same subclass over and over.
        After all, there is no difference if I put results to the namedtuple instance
        inside this function or outside the function by returning them.
        To make clear what is being stored, this function returns the result and let
        the outer function do the assignment.
    """

    if tmax is None:
        tmax = gg.star["time"].max()
    h, _ = np.histogram(gg.star["time"],weights=gg.star["m"]/(dt*1e9),bins=np.arange(tmin,tmax,dt))

    return h


def get_sfr(gg, dts, age_hist=None, dt_age_hist=None):
    """
        Or I can use age_hist if dts_this/dt_age_hist if integer.
        -> sum(age_hist[:dt_this/dt_age_hist]
        If not, this method will give an inaccurate answer.
    """
    if age_hist is None or dt_age_hist is None:
        try:
            if len(dts) > 0:
                sfrs= []
                for dt in dts:
                    sfrs.append(np.sum(gg.star["m"][gg.star["time"] < dt]))
                return sfrs
        except:
            # probably scalar.
            return np.sum(gg.star["m"][gg.star["time"] < dts])
    else:
        try:
            if len(dts) > 0:
                sfrs= []
                for dt in dts:
                    sfrs.append(np.sum(age_hist[:int(dt/dt_age_hist)]))

                return sfrs
        except:
            # probably scalar.
            return np.sum(age_hist[:int(dt/dt_age_hist)])


def get_galaxy_area(gg):
    """
        Returns the area coverd by stellar particles based on the gal.mmap.
    """
    #gg.npix_per_reff = 5
    dx = gg.rgal / gg.npix_per_reff # right?
    return dx*dx*sum(gg.mmap.ravel() > 0)
