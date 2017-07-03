import numpy as np

# 1. Gas properties
def get_gas_all(gg, CT, info, dr=5, rmax=200, density_ratio=1e-3):
    get_cold_cell(gg, info, dr, rmax, density_ratio)
    #mgas_tot, mgas_cold, Ln_gas = get_gas_properties(gg, info)
    gg.gas_results = CT.gas_results(get_gas_properties(gg, info))

def get_cold_cell(gg, info, dr=5, rmax=200, density_ratio=1e-3):
    """
        returns True if the galaxy has cold gas.
    """
    cold_cell = gg.cell[rho_t_cut(gg.cell, info)]
    if len(cold_cell) < 1:
        # this galaxy has no cold gas. That's possible. 
        return False
    else:
        #print("Numbe of cold cells", len(cold_cell))
        i_dense = ind_dense(cold_cell, dr=dr, rmax=rmax)
        if len(i_dense)/len(cold_cell) < density_ratio:
            print("Warning.. only a tiny fraction of cold gas is bound to the galaxy?")
        gg.cell = cold_cell[i_dense]
        return True
        #make_gal.extract_cold_gas(gg, dr=5, rmax=200)


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
    mgas_tot = np.sum(gg.cell[vrho]*gg.cell[vdx]**3)

    # Leave only cold gas and measure properties.
    has_cold_gas = get_cold_cell(gg, info)
    if has_cold_gas:
        vec_rot = np.cross(np.stack((gg.cell["x"],gg.cell["y"],gg.cell["z"])).T,
                           np.stack((gg.cell[vvx],
                                     gg.cell[vvy],
                                     gg.cell[vvz])).T)

        Ln_gas = (vec_rot.T * (gg.cell[vdx]**3 *gg.cell[vrho])).sum(axis=1)
    else:
        mgas_cold = 0
        Ln_gas = (-1,-1,-1)

    return mgas_tot, mgas_cold, Ln_gas

# 2. SFR
def get_sfr_all(gg, CT, sfr_dts=0.1,
                hist_dt=0.1, hist_tmin=0, hist_tmax=None,
                sfr_has_hist=False):
    h = get_age_hist(gg, dt=hist_dt, tmin=hist_tmin, tmax=hist_tmax)
    sfrs = get_sfr(gg, sfr_dts, has_hist=sfr_has_hist)
    area = get_galaxy_area(gg)

    gg.sfr_results= CT.sfr_results(dt, hist_tmin, hist_tmax, h/area, sfr_dts, sfrs/area, area)

def get_age_hiist(gg, 
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

    kpc
    if tmax is None:
        tmax = gg.star["time"].max()
    h, _ = np.histogram(gg.star["time"],weights=gg.star["m"]/(dt*1e9),bins=np.arange(tmin,tmax,dt))
    
    return h


def get_sfr(gg, dts, has_hist=False):
    """
        Or I can use age_hist if dts_this/dt_age_hist if integer.
        -> sum(age_hist[:dt_this/dt_age_hist]
        If not, this method will give an inaccurate answer.
    """
    if has_hist:
        return np.sum(age_hist[:int(dt_this/dt_age_hist)]
    else:
        try:
            if len(dts) > 0:
            sfrs= []
            for dt in dts:
                sfrs.append(np.sum(gg.star["m"][gg.star["time"] < dt]))
            return
        except:
            # probably scalar.
            return np.sum(gg.star["m"][gg.star["time"] < dts])


def get_galaxy_area(gg):
    """
        Returns the area coverd by stellar particles based on the gal.mmap.
    """
    dx = gg.rgal / gg.npix_per_reff # right?
    return dx*dx*sum(gg.mmap > 0)


##################

def get_vmax(gal):
    if gal.vmap is None:
        get_vmap(gal)
