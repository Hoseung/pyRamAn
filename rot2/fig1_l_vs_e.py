## Fig1.
import numpy as np
from numpy.lib.recfunctions import append_fields
from utils import match as mtc
import matplotlib.pyplot as plt

def density_map(x, y, sort=True):
    from scipy.stats import gaussian_kde
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    z /= max(z)

    idx = z.argsort()
    xx, yy = x[idx], y[idx]
    z = z[idx]

    #im = ax.scatter(xx, yy, c=z, s=50, edgecolor='')
    return xx,yy,z


def do_plot(x,y, obsdata,
            do_scatter=True,
            contour_label=False,
            surf = False,
            img_scale = 1.0,
            twocolors=['#4c72b0', '#c44e52'],
            den_cmap = "YlGnBu_r",
            levels=None,
            fname_vs_e = "./figs/lambda_vs_e_z0",
            d_alpha=1.0,
            sizeOfFont=12,
            label="This work",
            label2="",
            title="",
            ):
    #import scipy.stats as st

    fontsize_ticks = 6 * img_scale
    fontsize_tick_label = 8 * img_scale
    fontsize_legend = 5 * img_scale
    img_size_single_column =2.25 * img_scale

    #from matplotlib import rc, font_manager
#    fontProperties = {'family':'Liberation Sans',
#                      'weight' : 'normal', 'size' : sizeOfFont}
#    ticks_font = font_manager.FontProperties(family='Liberation Sans', style='normal',
#                   size=sizeOfFont, weight='normal', stretch='normal')
#    rc('text', usetex=True)
#    rc('font',**fontProperties)
    ticks_font = None

    xmin = ymin = -0.05
    xmax = ymax = 0.9

    fig, axmain=plt.subplots(1)

    fig.set_size_inches(img_size_single_column,
                        img_size_single_column)

    axmain.set_xlim(xmin, xmax)
    axmain.set_ylim(ymin, ymax)
    # suppress last tick
    axmain.set_xticks(np.arange(0, xmax, 0.1))
    axmain.set_yticks(np.arange(0, ymax, 0.1))
    axmain.set_xlabel(r"$\epsilon_{R_{eff}}$", fontsize=fontsize_tick_label, family="Liberation Sans")
    axmain.set_ylabel(r"$\lambda_{R_{eff}}$", fontsize=fontsize_tick_label)#, family="Liberation Sans")
    axmain.tick_params(axis='both', which='major', labelsize=fontsize_ticks)


    # S/R demarcation line
    sr_line_xx = np.arange(90)*0.01
    # Simple demarkation (Emsellem 2011)
    sr_line_yy = 0.31 * np.sqrt(sr_line_xx)
    axmain.plot(sr_line_xx, sr_line_yy, '--', lw=1, color='black')
    N_sr_em11 = sum(y < 0.31 * np.sqrt(x))

    # Cappellari+2016
    sr_line2_xx = np.arange(40)*0.01
    sr_line2_yy = 0.08+0.25*sr_line2_xx
    axmain.plot(sr_line2_xx, sr_line2_yy, lw=2, color='black')
    axmain.plot([0.4,0.4], [0, 0.08+0.1], lw=2, color='black')
    N_sr_ca16 = sum((y < 0.08+0.25*x) * (x < 0.4))

    Ngal = len(x)
    print("Number of galaxies in total", Ngal)
    print("SR(Emsellem 2011):", N_sr_em11)
    print("SR(Cappellari2016):", N_sr_ca16)
    axmain.text(0.6,0.1,
                "# FR: {} ({})".format(Ngal-N_sr_em11, Ngal-N_sr_ca16),
                fontsize=8)
    axmain.text(0.6,0.05,
                "# SR: {} ({})".format(N_sr_em11, N_sr_ca16),
                fontsize=8)

    # Draw main density map
    if surf:
        fname_vs_e = fname_vs_e + "_srf"
    elif contour_label:
        fname_vs_e = fname_vs_e + "_cl"
    try:
        nlevels = str(len(levels))
    except:
        nlevels = "default"
    fname_vs_e = fname_vs_e + "_" + nlevels
    xx,yy,z = density_map(x, y)
    axmain.scatter(xx, yy, c=z, s=15, edgecolor='',
                   cmap=den_cmap, rasterized=False,
                   alpha=1.0, label=label)
    axmain.set_xlim([-0.01,0.91])
    axmain.set_ylim([-0.01,0.91])

    if do_scatter:
        scatter = axmain.scatter(x,y, s=7,
                             facecolor=twocolors[0],
                             edgecolor='none',
                             alpha= 0.7,
                             label=label)
        fname_vs_e = fname_vs_e + "_sct"

    # My data
    if contour_label:
        axmain.clabel(cfset, inline=1, fontsize=7)

    #Observation data
    if obsdata is not None:
        axmain.scatter(obsdata[:,0], obsdata[:,1],
                   s=40,
                   color=twocolors[1],
                   marker=".",
                   lw=1,
                   alpha=0.8,
                   edgecolor='none',
                   label=label2)

    # Legend
    if 1 == 2:
        handles, labels = axmain.get_legend_handles_labels()
        #Create custom artists
        thisArtist = plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='')
        #Create legend from custom artist/label lists
        handles.append(thisArtist)
        labels.append(label)
        axmain.legend(handles, labels,
                      loc=2,
                      borderaxespad=0.,
                      labelspacing=1.2,
                      fontsize=fontsize_legend)
    else:
        axmain.legend(loc=2,
                      borderaxespad=0.,
                      labelspacing=1.2,
                      fontsize=fontsize_legend)

    plt.savefig(fname_vs_e + ".pdf", bbox_inches='tight')
    plt.savefig(fname_vs_e + ".png", bbox_inches='tight', dpi=200)
    axmain.set_title(title)
    plt.close()


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

    sd = sd[mtc.match_list_ind(common_ids, sd['id'])]
    dd = dd[mtc.match_list_ind(common_ids, dd['id'])]

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
            fname_vs_e = out_base+ "figs/"+prefix+"lambda_vs_e")
