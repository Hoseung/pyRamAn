
# coding: utf-8

# In[ ]:
import matplotlib.pyplot as plt
import numpy as np


# 1. Major Minor Rest contribution
def kde_den(data, cov=0.25):
    from scipy.stats import gaussian_kde
    density = gaussian_kde(data)
    xs = np.linspace(0,8,200)
    density.covariance_factor = lambda : cov
    density._compute_covariance()
    return density

def draw_kdes(dlM, dlm, dlo, dlt, ax
                  , nevents
                  , lw=1.5
                  , excess=False):

        dM = kde_den(dlM)
        dm = kde_den(dlm)
        do = kde_den(dlo)
        dtot = kde_den(dlt)

        xs=np.linspace(-0.7,0.7,51)
        if excess:
            i_positive = np.linspace(0.0,0.6,26)
            dM_curve = dM(xs) - np.concatenate((dM(i_positive)[::-1], dM(i_positive[1:])))
            dm_curve = dm(xs) - np.concatenate((dm(i_positive)[::-1], dm(i_positive[1:])))
            do_curve = do(xs) - np.concatenate((do(i_positive)[::-1], do(i_positive[1:])))
            dtot_curve = dtot(xs) - np.concatenate((dtot(i_positive)[::-1], dtot(i_positive[1:])))
        else:
            dM_curve = dM(xs)
            dm_curve = dm(xs)
            do_curve = do(xs)
            dtot_curve = dtot(xs)

        nM = len(dlM)
        nm = len(dlm)
        no = len(dlo)
        ntot = len(dlt)


        Mlabel="Major \n" +r"$N_{g}(N_{e})$" + " = {}({})".format(nM, nevents[0])
        mlabel="Minor \n" +r"$N_{g}(N_{e})$" + " = {}({})".format(nm, nevents[1])
        olabel="Rest  \n" +r"$N_{g}$" + " = {}".format(no)
        totlabel="Total\n" +r"$N_{g}$" + " = {}".format(ntot)
        ax.plot(xs, dM_curve*nM/ntot, label=Mlabel, lw=lw, color="r")
        ax.plot(xs, dm_curve*nm/ntot, label=mlabel, lw=lw, color="g")
        ax.plot(xs, do_curve*no/ntot, label=olabel, lw=lw, color="b")
        ax.plot(xs, dtot_curve, label=totlabel, lw=lw, color="black")

        ax.set_ylim([0, 1.15*ax.get_ylim()[1]])


def kde_sci(mpgs
    ,mstar_cut_hard = 5e9
    ,mcut=1e10
    ,fname="figs/test"
    ,wdir='./'
    ,nbins=21
    ,kde=True
    ,hist=True
    ,shade=True
    ,norm_hist=False
    ,pallette="muted"
    ,ylim=None
    ,per_event=True
    ,per_galaxy=True
    ,detected=True
    ,maj_ratio = 4
    ,excess=True
    ,img_scale=1.0
    ,examplegal=None):


    fontsize_ticks = 6 * img_scale
    fontsize_name = 8  * img_scale
    fontsize_tick_label = 8 * img_scale
    fontsize_legend = 5 * img_scale

    from matplotlib.ticker import NullFormatter

    l_dl_e = []
    l_mr_e = []
    l_mass_e = []

    s_dl_e = []
    s_mr_e = []
    s_mass_e = []

    l_dlt_g=[]
    l_dlo_g=[]
    l_dlM_g=[]
    l_dlm_g=[]
    l_mass_g=[]

    s_dlt_g=[]
    s_dlo_g=[]
    s_dlM_g=[]
    s_dlm_g=[]
    s_mass_g=[]


    M_changed = 0
    m_changed = 0
    no_merger_count = 0
    count = 0
    Maj_small = 0
    for i, gal in enumerate(mpgs):
        mgal = gal.data["mstar"][0]
        if mgal > mstar_cut_hard:
            delta_lambda_tot = np.average(gal.data['lambda_r'][:5]) - np.average(gal.data['lambda_r'][-5:])
            delta_lambda_major = 0
            delta_lambda_minor = 0

            # Large
            if mgal > mcut:
                if hasattr(gal, "merger"):
                    if gal.merger is not None:
                        l_dl_e.extend(gal.merger.delta_l)
                        l_mr_e.extend(gal.merger.mr)
                        for dl, mr in zip(gal.merger.delta_l, gal.merger.mr):
                            if (mr < maj_ratio) and (dl > -1):
                                delta_lambda_major = delta_lambda_major + dl
                            if (mr > maj_ratio) and (dl > -1):
                                delta_lambda_minor = delta_lambda_minor + dl

                delta_lambda_other = delta_lambda_tot - delta_lambda_major - delta_lambda_minor
                l_dlt_g.append(delta_lambda_tot)
                l_dlo_g.append(delta_lambda_other)
                l_dlM_g.append(delta_lambda_major)
                l_dlm_g.append(delta_lambda_minor)
            # small
            else:
                #s_mass_g.append(mgal)
                if hasattr(gal, "merger"):
                    if gal.merger is not None:
                        s_dl_e.extend(gal.merger.delta_l)
                        s_mr_e.extend(gal.merger.mr)
                        for dl, mr in zip(gal.merger.delta_l, gal.merger.mr):
                            if (mr < maj_ratio) and (dl > -1):
                                delta_lambda_major = delta_lambda_major + dl
                            if (mr > maj_ratio) and (dl > -1):
                                delta_lambda_minor = delta_lambda_minor + dl

                    delta_lambda_other = delta_lambda_tot - delta_lambda_major - delta_lambda_minor
                    s_dlt_g.append(delta_lambda_tot)
                    s_dlo_g.append(delta_lambda_other)
                    s_dlM_g.append(delta_lambda_major)
                    s_dlm_g.append(delta_lambda_minor)

    l_dlt_g = np.array(l_dlt_g)
    l_dlo_g = np.array(l_dlo_g)
    l_dlM_g = np.array(l_dlM_g)
    l_dlm_g = np.array(l_dlm_g)
    #l_mass_g = np.array(l_mass_g)

    s_dlt_g = np.array(s_dlt_g)
    s_dlo_g = np.array(s_dlo_g)
    s_dlM_g = np.array(s_dlM_g)
    s_dlm_g = np.array(s_dlm_g)
    #s_mass_g = np.array(s_mass_g)

    # detected
    l_dlM_g = l_dlM_g [l_dlM_g !=0]
    #l_dlM_M = l_mass_g[l_dlM_g !=0]
    l_dlm_g = l_dlm_g [l_dlm_g !=0]
    #l_dlm_M = l_mass_g[l_dlm_g !=0]
    #l_dlo_M = l_mass_g

    s_dlM_g = s_dlM_g [s_dlM_g !=0]
    #s_dlM_M = s_mass_g[s_dlM_g !=0]
    s_dlm_g = s_dlm_g [s_dlm_g !=0]
    #s_dlm_M = s_mass_g[s_dlm_g !=0]
    #s_dlo_M = s_mass_g


    l_dl_e = np.array(l_dl_e)
    l_mr_e = np.array(l_mr_e)
    #l_mass_e = []

    s_dl_e = np.array(s_dl_e)
    s_mr_e = np.array(s_mr_e)
    #s_mass_e = []

    fig, axs = plt.subplots(3, sharex=True)
    fig.set_size_inches(4.75,7)
    plt.subplots_adjust(hspace=0.01)

    all_dlM_g = np.concatenate((l_dlM_g,s_dlM_g))
    all_dlm_g = np.concatenate((l_dlm_g,s_dlm_g))
    all_dlo_g = np.concatenate((l_dlo_g,s_dlo_g))
    all_dlt_g = np.concatenate((l_dlt_g,s_dlt_g))

    draw_kdes(all_dlM_g,
              all_dlm_g,
              all_dlo_g,
              all_dlt_g,
              axs[0],
              [sum(s_mr_e < maj_ratio) + sum(l_mr_e < maj_ratio),
               sum(s_mr_e > maj_ratio) + sum(l_mr_e > maj_ratio),
               len(all_dlo_g)],
              excess=excess)

    draw_kdes(l_dlM_g,
              l_dlm_g,
              l_dlo_g,
              l_dlt_g,
              axs[1],
              [sum(l_mr_e < maj_ratio),
               sum(l_mr_e > maj_ratio),
               len(l_dlo_g)],
              excess=excess)

    draw_kdes(s_dlM_g,
              s_dlm_g,
              s_dlo_g,
              s_dlt_g,
              axs[2],
              [sum(s_mr_e < maj_ratio),
               sum(s_mr_e > maj_ratio),
               len(s_dlo_g)],
              excess=excess)

    axs[0].set_xlim([-0.6,0.6])
    for ax in axs:
        ax.xaxis.grid()
        leg = ax.legend(fontsize=fontsize_legend)
        leg.get_frame().set_alpha(0.5)
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.set_ylabel("relative probability", fontsize=fontsize_tick_label)

    axs[2].set_xlabel(r"$\Delta \lambda_{R_{eff}}$", fontsize=fontsize_tick_label, family="Liberation Sans")
    axs[2].tick_params(labelsize=fontsize_ticks)
    axs[2].set_xlim([-0.7,0.6])
    #axs[0].legend(fontsize=12)

    axs[0].text(0.05, 0.87, "(A)", weight="bold", transform=axs[0].transAxes, fontsize=fontsize_name)
    axs[0].text(0.15, 0.87, "All",transform=axs[0].transAxes, fontsize=fontsize_name)
    axs[1].text(0.05, 0.87, "(B) ", weight="bold",transform=axs[1].transAxes, fontsize=fontsize_name)
    axs[1].text(0.15, 0.87, r"$log_{10}M_{\star} > $ " +"{:.1f}".format(np.log10(mcut))
                , fontsize=fontsize_name
                , transform=axs[1].transAxes)
    axs[2].text(0.05, 0.87, "(C) ", weight="bold",transform=axs[2].transAxes, fontsize=fontsize_name)
    axs[2].text(0.15, 0.87, r"$log_{10}M_{\star} < $ " +"{:.1f}".format(np.log10(mcut))
                , fontsize=fontsize_name
                , transform=axs[2].transAxes)

    if examplegal is not None:
        dls = examplegal.merger.delta_l
        mrs = examplegal.merger.mr

        for dl, mr in zip(dls, mrs):
            if mr < maj_ratio:
                axs[0].scatter(dl, [0.15], facecolor='r', edgecolor="w", marker="d", s=40, zorder=20)
            else:
                axs[0].scatter(dl, [0.15], facecolor='g', edgecolor="w", marker="d", s=40, zorder=20)

        dl_tot = examplegal.data["lambda_r"][0] - examplegal.data["lambda_r"][-1]
        axs[0].scatter(dl_tot, [0.15], facecolor='black', edgecolor="w", marker="d", s=40, zorder=20)
        dl_o = dl_tot - sum(examplegal.merger.delta_l)
        axs[0].scatter(dl_o, 0.15, facecolor='b', edgecolor="w", marker="d", s=40, zorder=20)

        # legend
        axs[0].scatter(0.21, 1, facecolor='none', edgecolor="black", marker="d", s=40)
        axs[0].text(0.24, 0.95, "example galaxy", fontsize=8)
        print("Example gal")


    plt.savefig(fname + "{:.1f}.png".format(np.log10(mcut)), dpi=200, bbox_inches="tight")
    plt.savefig(fname + "{:.1f}.pdf".format(np.log10(mcut)), bbox_inches='tight') # eps does NOT support transparency!
    plt.savefig(fname + "{:.1f}.eps".format(np.log10(mcut)), bbox_inches='tight')
    plt.savefig(fname + "{:.1f}.svg".format(np.log10(mcut)), bbox_inches='tight')

    plt.close()



# In[ ]:



# 2. Lambda Vs. Ellipticity
def plot_density_map(axmain, x, y, xmin, xmax, ymin, ymax,
                    levels=None, color=True, cmap="winter",
                    surf=False, bw_method="silverman",
                    d_alpha=1.0):
    import scipy.stats as st
    # Draw main density map
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])

    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values, bw_method=bw_method)
    f = np.reshape(kernel(positions).T, xx.shape)
    # , [0.2, 2, 3, 5, 10, 16]  custom contour levels.
    f /= max(f.ravel())

#   ains = inset_axes(axmain, width='5%', height='60%', loc=5)
    if surf:
        cfset = axmain.contourf(xx, yy, f,
                            levels=levels,
                            cmap=cmap,
                            alpha = d_alpha)#,

    else:
        cfset = axmain.contour(xx, yy, f,
                            levels=levels,
                            cmap=cmap,
                            alpha = d_alpha,
                            linewidths=0.6)

    return cfset
def plot_sami(ax, data, contour=True, scatter=False):
    x = data['ellp']
    y = data['r1']

    xmin, xmax = -0.05, 0.8
    ymin, ymax = -0.05, 0.8

    if contour:
        plot_densiy_map(ax, x, y, xmin, xmax, ymin, ymax,
                        N=4,
                        levels=None, color=False, cmap="spring")
    if scatter:
        ax.scatter()

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

def do_plot(x,y, atlas,
            do_scatter=True,
            contour_label=False,
            surf = False,
            img_scale = 1.0,
            twocolors=['#4c72b0', '#c44e52'],
            den_cmap = "PuBu",
            levels=None,
            fname_vs_e = "./figs/lambda_vs_e_z0",
            d_alpha=1.0,
            sizeOfFont=12
            ):
    import scipy.stats as st

    fontsize_ticks = 6 * img_scale
    fontsize_tick_label = 8 * img_scale
    fontsize_legend = 5 * img_scale
    img_size_single_column =2.25 * img_scale

    from matplotlib import rc, font_manager
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
#    for label in axmain.get_xticklabels():
#        label.set_fontproperties(ticks_font)
#    for label in axmain.get_yticklabels():
#        label.set_fontproperties(ticks_font)
    #axmain.set_xticklabels(axmain.get_xticks(), family="helvetica")
    #axmain.set_yticklabels(axmain.get_yticks(), family="helvetica")

    # S/R demarcation line
    sr_line_xx = np.arange(90)*0.01
    # Simple demarkation (Emsellem 2011)
    sr_line_yy = 0.31 * np.sqrt(sr_line_xx)
    axmain.plot(sr_line_xx, sr_line_yy, '--', lw=1, color='black')

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

    if 1 == 2:
        cfset = plot_density_map(axmain, x, y, xmin, xmax, ymin, ymax,
                            levels=levels,
                            cmap=den_cmap,
                            surf=True,
                            d_alpha=d_alpha)

    if 1==1:
        xx,yy,z = density_map(x, y)
        axmain.scatter(xx, yy, c=z, s=15, edgecolor='',
                       cmap=den_cmap, rasterized=False,
                       alpha=1.0, label="This work")

    if do_scatter:
        scatter = axmain.scatter(x,y, s=7,
                             facecolor=twocolors[0],
                             edgecolor='none',
                             alpha= 0.7,
                             label="This work")
        fname_vs_e = fname_vs_e + "_sct"

    if 1 == 3:
        cfset = plot_density_map(axmain, x, y, xmin, xmax, ymin, ymax,
                            levels=levels,
                            cmap="winter",
                            surf=False,
                            d_alpha=d_alpha)
    # My data


    if contour_label:
        axmain.clabel(cfset, inline=1, fontsize=7)

    #ATLAS3D
    axmain.scatter(atlas[:,0], atlas[:,1],
                   s=40,
                   color=twocolors[1],
                   marker=".",
                   lw=1,
                   alpha=0.8,
                   edgecolor='none',
                   label="ATLAS" + r"$^{3D}$")

    # Legend
    if 1 == 2:
        handles, labels = axmain.get_legend_handles_labels()
        #Create custom artists
        thisArtist = plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='')
        #Create legend from custom artist/label lists
        handles.append(thisArtist)
        labels.append("This work")
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
    #plt.savefig(fname_vs_e + ".svg", bbox_inches='tight')
    #plt.savefig(fname_vs_e + ".eps", bbox_inches='tight')

    #plt.show()

    plt.close()

def truncate_colormap(cmap_name, minval=0.0, maxval=1.0, n=100):
    """
    Use only part of color maps
    """
    import matplotlib.colors as colors
    cmap = plt.get_cmap(cmap_name)
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap




# 3. Lambda evolution


from astropy.cosmology import WMAP7, z_at_value

class Merger():
    pass

def aexp2zred(aexp):
    return [1.0/a - 1.0 for a in aexp]

def zred2aexp(zred):
    return [1.0/(1.0 + z) for z in zred]

def lbt2aexp(lts):
    import astropy.units as u
    zreds = [z_at_value(WMAP7.lookback_time, ll * u.Gyr) for ll in lts]
    return [1.0/(1+z) for z in zreds]


def age2zred(lts):
    import astropy.units as u
    zreds = [z_at_value(WMAP7.age, ll * u.Gyr) for ll in lts]
    return zreds

def l_at_smoothed_r(gal, npix_per_reff=5):
    import analysis.Major_Minor_accretion as mma
    n_valid_points = sum(gal.data["reff"] > 0)
    new_l_arr = np.zeros(n_valid_points)
    new_reff = mma.smooth(gal.data["reff"])#[gal.data["reff"] > 0])
    for i in range(n_valid_points):
        try:
            lambdar = gal.data["lambda_arr"][i]
            ind_org = npix_per_reff - 1
            #i_new =  # 0-indexed.
            ind_new = new_reff[i]/gal.data["reff"][i] * ind_org
            il = np.fix(ind_new).astype(int)
            ir = il + 1
            if ir >= len(lambdar):
                new_l_arr[i] = lambdar[-1]
            else:
                new_l_arr[i] = lambdar[il]*(ir-ind_new) + lambdar[ir]*(ind_new-il)
                #new_l_arr[i] = lambdar[il]*(ir-ind_org) + lambdar[ir]*(ind_org-il)
        except:
            new_l_arr[i] = gal.data["lambda_arr"][i] # = 0 with bad measurements.
    return new_l_arr


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


def modify_ticks(zreds, aexps, ax, ax2, nout_ini, nout_fi):
    modify_ticks1(zreds, aexps, ax, nout_ini, nout_fi)
    modify_ticks2(zreds, aexps, ax2, nout_ini, nout_fi)

def modify_ticks1(zreds, aexps, ax, nout_ini, nout_fi, fontsize=12):
    nbins = 40 # along the y-axis.
    nnouts = nout_fi - nout_ini + 1

    # For a given list of nouts,
    # calculate a nice-looking set of zreds AND lookback times
    zz_target = [0, 0.5, 1.0, 2.0, 3.0]#[::-1]
    x_tick_pos =np.array([187, 120, 87, 54, 37]) - nout_ini
    #x_tick_pos = np.searchsorted(zreds[::-1], zz_target)[::-1]# + nout_ini# + nout_min
    # searchsorted requires arrays be in ascending order.

    ax.set_xlabel("Redshift", fontsize=fontsize, family="Liberation Sans")
    ax.set_xlim([0, nnouts])

    ax.set_xticks(x_tick_pos)#[::-1])
    ax.set_xticklabels(labels = ["{:0.1f}".format(z) for z in zz_target])


def modify_ticks2(zreds, aexps, ax2, nout_ini, nout_fi):
    nbins = 40 # along the y-axis.
    nnouts = nout_fi - nout_ini + 1

    # For a given list of nouts,
    # calculate a nice-looking set of zreds AND lookback times
    zz_target = [0, 0.2, 0.5, 1.0, 2.0, 3.0]#[::-1]
    x_tick_pos = np.searchsorted(zreds[::-1], zz_target)[::-1]# + nout_ini# + nout_min
    # searchsorted requires arrays be in ascending order.

    #  at z = 3, lbt = 11.5243, age of universe = 2.18844
    u_age_targets=[3,5,8,10,13]
    u_age_target_str=["{:.0f}".format(l) for l in u_age_targets]
    z_targets_u_age = age2zred(u_age_targets)
    u_age_pos = np.searchsorted(zreds[::-1], z_targets_u_age)[::-1]

    ax2.set_xticks(u_age_pos)
    ax2.set_xticklabels(labels = u_age_target_str)
    ax2.set_xlabel("Age of the universe (Gyr)", family="Liberation Sans")




def plot_lambda_evol3(mpgs, fig, axs, nout_ini, nout_fi,
                     wdir_info='./',
                     density = "none",
                     cmap ="jet",
                     img_scale=1.0,
                     sizeOfFont=9):

    fontsize_ticks = 6 * img_scale
    fontsize_tick_label = 8 * img_scale
    fontsize_legend = 5 * img_scale

    from matplotlib import rc, font_manager
    fontProperties = {'family':'Liberation Sans',
                      'weight' : 'normal', 'size' : sizeOfFont}
    ticks_font = font_manager.FontProperties(family='Liberation Sans', style='normal',
                   size=sizeOfFont, weight='normal', stretch='normal')
    rc('text', usetex=True)
    rc('font',**fontProperties)


    nnouts = nout_fi - nout_ini + 1
    ngals_tot = len(mpgs)
    lambda_evol_all = np.zeros([ngals_tot, nnouts])

    # Starting nouts of galaxies are different.
    # But all end at nout = 187.
    for igal, gal in enumerate(mpgs):
        for inout, nout in enumerate(gal.nouts):
            lambda_evol_all[igal][nout - nout_ini] = gal.data['lambda_r'][inout]


    # change ticks
    zreds=[]
    aexps=[]

    # Only nout_ini < nout < nout_fi values are taken.
    # So are the x ticks.
    import load
    for nout in range(nout_ini, nout_fi + 1):
        info = load.info.Info(nout=nout, base=wdir_info, load=True)
        aexps.append(info.aexp)
        zreds.append(info.zred)
    aexps = np.array(aexps)
    zreds = np.array(zreds)

    modify_ticks1(zreds, aexps, axs[2], nout_ini, nout_fi, fontsize=fontsize_tick_label)

    if density == "hexbin":
        xx = np.tile(np.arange(nnouts), ngals_tot)
        all_data = lambda_evol_all.ravel()
        ind_ok = all_data > 0.01
        #xx,yy,z = density_map(xx[ind_ok], all_data[ind_ok])
        for ax in axs:
            im = ax.hexbin(xx[ind_ok], all_data[ind_ok],
                       gridsize=int(nnouts/3),
                       bins=None,
                       cmap=cmap)
            yticks_ok=[0.0, 0.2, 0.4, 0.6, 0.8]
            ax.set_ylim([-0.1, 0.9])
            ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
            ax.set_yticklabels([str(yy) for yy in yticks_ok])
            ax.set_ylabel(r"$\lambda_{R_{eff}}$", fontsize=fontsize_tick_label, family="Liberation Sans")
            ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks, family="Liberation Sans")

    elif density == "kernel":
        xx = np.tile(np.arange(nnouts), ngals_tot)
        all_data = lambda_evol_all.ravel()

        ind_ok = all_data > 0.01
        lambda_range=[0.01, 0.8]
        xx,yy,z = density_map(xx[ind_ok], all_data[ind_ok])

        for ax in axs:
            ax.scatter(xx, yy, c=z, s=50, edgecolor='', cmap=cmap, rasterized=True)
#            yticks_ok=[0.0, 0.2, 0.4, 0.6, 0.8]
#            ax.set_ylim([-0.05, 0.9])
#            ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
#            ax.set_yticklabels([str(yy) for yy in yticks_ok])
#            ax.set_ylabel(r"$\lambda_{R_{eff}}$", fontsize=fontsize_tick_label, family="Liberation Sans")
#            ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)

#    axs[2].tick_params(axis='x', which='major', labelsize=fontsize_ticks)


    return zreds, aexps



def simple_evol_path(gg, dt=10):
    """
    return mass and lambda at point (by taking averages?)
    """
    nout = []
    mass = []
    lamb = []

    for i in range(len(gg.nouts)-1, 0, -dt):
        nout.append(gg.nouts[i])
        mass.append(gg.data["mstar"][i])
        lamb.append(gg.smoothed_lambda[i])
    return np.array(nout), np.array(mass), np.array(lamb)

def plot_minor(twogals, ax, suptitle="",
                  normalize=True,
                  style="arrow",
                  nout_ini=37,
                  img_scale=2.0,
                  annotate="(B )",
                  arrow_scale = 40):
    import matplotlib as mpl
    """
    Up to 4 galaxies in a plot
    """
    #im = pickle.load(open("lambda_evol_all.pickle", "rb"))
    #im.axes.set_ybound([-0.1, 0.85])

    fontsize_ticks = 8 * img_scale
    fontsize_tick_label = 8 * img_scale
    fontsize_legend = 7 * img_scale
    img_size_single_column =2.25 * img_scale

    loc = [0, 0, 0, 0]
    if len(twogals) == 2:
        if twogals[0].smoothed_lambda[-1] < twogals[1].smoothed_lambda[-1]:
            loc[1] = 0.65
        else:
            loc[0] = 0.65

    #ax.set_title(suptitle)

    for igal, gg in enumerate(twogals):
        if len(gg.smoothed_lambda) > 100:
            nn, mm, ll = simple_evol_path(gg)
            nn -= nout_ini
            if normalize:
                widths = (mm[:-1]/mm[-1])
            else:
                widths = np.sqrt(mm[:-1]/1e9)
            for i in range(len(mm)-1):
                #mm = im.axes.quiver(nn[i], ll[i], nn[i+1]-nn[i], ll[i+1]-ll[i],
                if style=="arrow":
                    ax.arrow(nn[i], ll[i], nn[i+1]-nn[i], ll[i+1]-ll[i],
                             fc="none",
                             ec=["red", "green", "yellow", "white"][igal],
                             alpha=0.5,
                             width=widths[i],
                             head_width=10 * widths[i],
                             head_length=2000 * widths[i])

                elif style=="quiver":
                    #ax.quiver(nn[:-1], ll[:-1], nn[1:]-nn[:-1], ll[1:]-ll[:-1],
                    ax.quiver(nn[i], ll[i], nn[i+1]-nn[i], ll[i+1]-ll[i],
                           scale_units='xy', angles='xy', scale=1,
                           linewidth = 0.5,
                           headwidth = 2 * widths[i]**0.5,
                           headlength = 0.5 * widths[i],
                           #headaxislength = 0.3 * widths[i],
                           alpha=0.3, facecolor="none",
                           edgecolor=["red", "green", "yellow", "white"][igal])
                          #linewidth = 0.005 * widths[i]**3,
                    # distance between two points differ since x distance is fixed but
                    # y distance changes. So fixed headeraxislength results in changing
                    # head shape. It's not trivial to retain the shape of arrows.
                elif style=="stream":
                    patch = mpl.patches.FancyArrowPatch(
                        (nn[i], ll[i]), (nn[i+1],ll[i+1]),
                        arrowstyle=['->','-|>',"fancy","simple","wedge","-"][3],
                        mutation_scale=np.sqrt(widths[i]) * arrow_scale,
#                        lw=np.sqrt(widths[i]) * 5,
                        edgecolor=["red", "green", "yellow", "white"][igal],
                        facecolor="white",
                        fill=True,
                        joinstyle=['miter', 'round', 'bevel'][0],
                        alpha = 0.5,
                        arrow_transmuter=None,
                        connectionstyle='Arc3',
                        connector=None,
                        patchA=None,
                        patchB=None,
                        shrinkA=2.0,
                        shrinkB=2.0,
                        mutation_aspect=None,
                        dpi_cor=1.0)
                    ax.add_patch(patch)
                    # Draw edge once more to have vivid edge + transparent face
                    patch = mpl.patches.FancyArrowPatch(
                        (nn[i], ll[i]), (nn[i+1],ll[i+1]),
                        arrowstyle=['->','-|>',"fancy","simple","wedge","-"][3],
                        mutation_scale=np.sqrt(widths[i]) * arrow_scale,
                        edgecolor=["red", "green", "yellow", "white"][igal],
                        fill=False,
                        joinstyle=['miter', 'round', 'bevel'][0],
                        alpha = 0.3,
                        arrow_transmuter=None,
                        connectionstyle='Arc3',
                        connector=None,
                        patchA=None,
                        patchB=None,
                        shrinkA=2.0,
                        shrinkB=2.0,
                        mutation_aspect=None,
                        dpi_cor=1.0)
                    ax.add_patch(patch)

            if gg.merger is not None:
                for nouti, nout, mr in zip(gg.merger.nout_ini, gg.merger.nout, gg.merger.mr):
                    #im.axes.axvline(nouti - nout_ini, linestyle=':')
                    #im.axes.axvline(nout - nout_ini, linestyle=':')
                    ypos_annotate = gg.smoothed_lambda[gg.nouts == nout] + 0.05
                    ax.plot([nout-nout_ini, nout-nout_ini],
                                 [loc[igal], ypos_annotate - 0.1 + loc[igal] * 0.1/0.75],
                                 lw = 1,
                                 linestyle="--",
                                 color = ["red", "green", "yellow", "white"][igal],
                                 alpha = 0.5)
                    ax.annotate("{:.2f}".format(1./mr),
                                    xy=(nout - nout_ini, 0.32 + ((loc[igal]) - 0.32) * 1.1),
                                    fontsize=fontsize_legend)
                    # text slightly further than the end of dashed lines
            #ax.annotate(suptitle, xy=(5, 0.79), fontsize=fontsize_ticks)
        else:
            print("Too short")
    ax.text(0.05, 0.9, annotate, weight="bold",transform=ax.transAxes, fontsize=fontsize_ticks)
    ax.text(0.13, 0.9, suptitle, transform=ax.transAxes, fontsize=fontsize_ticks)


def plot_rest(twogals, ax,
          suptitle="",
          normalize=True,
          style="arrow",
          nout_ini=37,
          img_scale=2.0,
          annotate="(B )",
          arrow_scale = 40):

    import matplotlib as mpl
    """
    Up to 4 galaxies in a plot
    """
    #im = pickle.load(open("lambda_evol_all.pickle", "rb"))
    #im.axes.set_ybound([-0.1, 0.85])

    fontsize_ticks = 8 * img_scale
    fontsize_tick_label = 8 * img_scale
    #fontsize_legend = 5 * img_scale
    img_size_single_column =2.25 * img_scale

    loc = [0, 0, 0, 0]
    if len(twogals) == 2:
        if twogals[0].smoothed_lambda[-1] < twogals[1].smoothed_lambda[-1]:
            loc[1] = 0.65
        else:
            loc[0] = 0.65

    #ax.set_title(suptitle)

    for igal, gg in enumerate(twogals):
        if len(gg.smoothed_lambda) > 100:
            nn, mm, ll = simple_evol_path(gg)
            nn -= nout_ini
            if normalize:
                widths = (mm[:-1]/mm[-1])
            else:
                widths = np.sqrt(mm[:-1]/1e9)
            for i in range(len(mm)-1):
                patch = mpl.patches.FancyArrowPatch(
                    (nn[i], ll[i]), (nn[i+1],ll[i+1]),
                    arrowstyle=['->','-|>',"fancy","simple","wedge","-"][3],
                    mutation_scale=np.sqrt(widths[i]) * arrow_scale,
#                    lw=np.sqrt(widths[i]) * 5,
                    edgecolor=["red", "green", "yellow", "white"][igal],
                    facecolor="white",
                    fill=True,
                    joinstyle=['miter', 'round', 'bevel'][0],
                    alpha = 0.5,
                    arrow_transmuter=None,
                    connectionstyle='Arc3',
                    connector=None,
                    patchA=None,
                    patchB=None,
                    shrinkA=2.0,
                    shrinkB=2.0,
                    mutation_aspect=None,
                    dpi_cor=1.0)
                ax.add_patch(patch)
                # Draw edge once more to have vivid edge + transparent face
                patch = mpl.patches.FancyArrowPatch(
                    (nn[i], ll[i]), (nn[i+1],ll[i+1]),
                    arrowstyle=['->','-|>',"fancy","simple","wedge","-"][3],
                    mutation_scale=np.sqrt(widths[i]) * arrow_scale,
                    edgecolor=["red", "green", "yellow", "white"][igal],
                    fill=False,
                    joinstyle=['miter', 'round', 'bevel'][0],
                    alpha = 0.3,
                    arrow_transmuter=None,
                    connectionstyle='Arc3',
                    connector=None,
                    patchA=None,
                    patchB=None,
                    shrinkA=2.0,
                    shrinkB=2.0,
                    mutation_aspect=None,
                    dpi_cor=1.0)
                ax.add_patch(patch)
        else:
            print("Too short")
    ax.text(0.05, 0.9, annotate, weight="bold",transform=ax.transAxes, fontsize=fontsize_ticks)
    ax.text(0.13, 0.9, suptitle, transform=ax.transAxes, fontsize=fontsize_ticks)


def plot_major(twogals, ax,
                  suptitle="",
                  normalize=True,
                  nout_ini=37,
                  img_scale=2.0,
                  arrow_scale = 40):
    import matplotlib as mpl

    fontsize_ticks = 8 * img_scale
    fontsize_tick_label = 8 * img_scale
    fontsize_legend = 7 * img_scale

    loc = [0, 0.65, 0, 0]

    for igal, gg in enumerate(twogals):
        if len(gg.smoothed_lambda) > 100:
            nn, mm, ll = simple_evol_path(gg)
            nn -= nout_ini
            if normalize:
                widths = (mm[:-1]/mm[-1])
            else:
                widths = np.sqrt(mm[:-1]/1e9)
            for i in range(len(mm)-1):
                patch = mpl.patches.FancyArrowPatch(
                    (nn[i], ll[i]), (nn[i+1],ll[i+1]),
                    arrowstyle=['->','-|>',"fancy","simple","wedge","-"][3],
                    mutation_scale=np.sqrt(widths[i]) * arrow_scale,
#                        lw=np.sqrt(widths[i]) * 5,
                    edgecolor=["red", "green", "yellow", "white"][igal],
                    facecolor="white",
                    fill=True,
                    joinstyle=['miter', 'round', 'bevel'][0],
                    alpha = 0.5,
                    dpi_cor=1.0)
                ax.add_patch(patch)
                # Draw edge once more to have vivid edge + transparent face
                patch = mpl.patches.FancyArrowPatch(
                    (nn[i], ll[i]), (nn[i+1],ll[i+1]),
                    arrowstyle=['->','-|>',"fancy","simple","wedge","-"][3],
                    mutation_scale=np.sqrt(widths[i]) * arrow_scale,
                    edgecolor=["red", "green", "yellow", "white"][igal],
                    fill=False,
                    joinstyle=['miter', 'round', 'bevel'][0],
                    alpha = 0.3,
                    dpi_cor=1.0)
                ax.add_patch(patch)

            if gg.merger is not None:
                for nouti, nout, mr in zip(gg.merger.nout_ini, gg.merger.nout, gg.merger.mr):
                    if mr < 4:
                    #im.axes.axvline(nouti - nout_ini, linestyle=':')
                    #im.axes.axvline(nout - nout_ini, linestyle=':')
                        ypos_annotate = gg.smoothed_lambda[gg.nouts == nout] + 0.05
                        ax.plot([nout-nout_ini, nout-nout_ini],
                                     [loc[igal], ypos_annotate - 0.1 + loc[igal] * 0.1/0.75],
                                     lw = 1,
                                     linestyle="--",
                                     color = ["red", "green", "yellow", "white"][igal],
                                     alpha = 0.5)
                        ax.annotate("{:.2f}".format(1./mr),
                                    xy=(nout - nout_ini, 0.32 + ((loc[igal]) - 0.32) * 1.1),
                                    fontsize=fontsize_legend)
                        # text slightly further than the end of dashed lines

        else:
            print("Too short")
    ax.text(0.05, 0.9, "(A) ", weight="bold",transform=ax.transAxes, fontsize=fontsize_ticks)
    ax.text(0.13, 0.9, suptitle, transform=ax.transAxes, fontsize=fontsize_ticks)
