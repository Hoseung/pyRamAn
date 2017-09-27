## Fig1.
import numpy as np
from numpy.lib.recfunctions import append_fields
from analysis.all_plot_modules import do_plot
from utils import match as mtc

def lambda_vs_e(gal_data,
                data_type="allresults",
                istep=0,
                mass_bin=None,
                out_base='./',
                prefix=""):
    """
    Parameters
    ----------
        gal_data : allresults at a nout. (allresults[0])
        date_type : "allresults" or "serial_results"
        mass_bin : list   [1,1e12]
            If not None, only the galaxies inside the mass_bin are considered.
    """
    # OBS data points (Lambda and eps)
    obs_dir='/home/hoseung/Work/data/Rotator_data/'
    # load SAMI data

    sd = np.genfromtxt(obs_dir+'/samidata/jvds_sami_lambdar_220216_0arc_fwhm.txt',
                        dtype=(int, float, float, float, float, float, float),
                        usecols=(0,2,3,6,7,10,11))
    sd.dtype.names = ('id', 'r05', 'e05', 'r1', 'e1', 'r2', 'e2') # e = errors, ellip = ellipticity

    dd = np.genfromtxt(obs_dir+'/samidata/input_catalog_lambdar.txt',
                              dtype=(int, float, float, float),
                              usecols=(0,3,4,6))
    dd.dtype.names = ('id', 're', 'ellp', 'mstar')

    common_ids = np.intersect1d(sd['id'], dd['id'])

    ind = mtc.match_list_ind(common_ids, sd['id'])
    sd = sd[ind]
    ind = mtc.match_list_ind(common_ids, dd['id'])
    dd = dd[ind]

    if all(sd['id'] == dd['id']):
        tt = dd[['re', 'ellp', 'mstar']]
        sami_data = append_fields(sd, tt.dtype.names, [tt[n] for n in tt.dtype.names])
        # only valid values
        all_finite = np.all([np.isfinite(sami_data[field]) for field in ['id', 'r05', 'r1', 're', 'ellp', 'mstar']], axis=0)
        sami_data = sami_data[all_finite]

    # exclude ellp == 0
    sami_data = sami_data[sami_data['ellp'] > 0].data

    # Mass cut
    sami_data = sami_data[sami_data['mstar'] > 10.0]

    atlas = np.genfromtxt(obs_dir+'ATLAS3D/Emsellem2011_Atlas3D_Paper3_TableB1.txt',
                          skip_header=12,
                          usecols=(2,7))

    # Simulation data
    lambda_e=[]
    eps=[]
    mass=[]

    if data_type == "allresults":
        for gg in gal_data:
            try:
                lambda_e.append(gg.lambda_r[0])
                eps.append(gg.mge_result_list[0]["eps"])
                mass.append(gg.mstar)
            except:
                pass
    elif data_type == "serial_results":
        for gg in gal_data:
            try:
                lambda_e.append(gg.finedata["lambda_r"][istep])
                eps.append(gg.finedata["eps"][istep])
                mass.append(gg.finedata["mstar"][istep])
            except:
                pass
    elif data_type == "array":
        lambda_e = gal_data["lambda_r"][:,istep]
        eps = gal_data["eps"][:,istep]
        mass = gal_data["mstar"][:,istep]

    i = np.isfinite(lambda_e) * lambda_e > 0
    x = np.array(eps)[i] # isfinit
    y = np.array(lambda_e)[i]
    mass = np.array(mass)[i]

    if mass_bin is not None:
        x = x[(mass > mass_bin[0]) * ( mass < mass_bin[1])]
        y = y[(mass > mass_bin[0]) * ( mass < mass_bin[1])]
        mass = mass[(mass > mass_bin[0]) * ( mass < mass_bin[1])]

    print("number of galaxies in total", len(x))
    print("Number of galaxies below the demarkation line:", sum(y < 0.31 * np.sqrt(x)))

    # Color
    twocolors=['#4c72b0', '#c44e52']
    #red = '#c44e52'
    do_plot(x,y, atlas,
            do_scatter=False,
            contour_label=False,
            surf = False,
            img_scale = 2.0,
            twocolors=twocolors,
            den_cmap = "YlGnBu_r",
            d_alpha=1.0,
            levels=None,#np.linspace(0.02, 1.0, 19),
            fname_vs_e = out_base+ "figs/"+prefix+"lambda_vs_e_HM")
    do_plot(x,y, -10*np.arange(6).reshape(2,3),
            do_scatter=False,
            contour_label=False,
            surf = False,
            img_scale = 2.0,
            twocolors=twocolors,
            den_cmap = "YlGnBu_r",
            d_alpha=1.0,
            levels=None,#np.linspace(0.02, 1.0, 19),
            fname_vs_e = out_base+ "figs/"+prefix+"lambda_vs_e_")
