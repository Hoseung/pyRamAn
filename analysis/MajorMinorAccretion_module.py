import matplotlib.pyplot as plt
import numpy as np
import analysis.evol_lambda as evl
import analysis.Major_Minor_accretion as mma
import analysis.misc as amsc
import tree.ctutils as ctu
import utils.match as mtc


def gaussian_fit(ax, data, dx, color='green'):
    import matplotlib.mlab as mlab
    mean = np.mean(data)
    variance = np.var(data)
    sigma = np.sqrt(variance)
    x = np.linspace(min(data), max(data), 100)
    scale = len(data)*dx
    ax.plot(mlab.normpdf(x,mean,sigma)*scale, x, color=color)
    ax.text(0.1, 0.9, "mean = {:.2f}\n sig = {:.2f}".format(mean,sigma),
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes)

def fwhm(xx, curve):
    from scipy.interpolate import UnivariateSpline
    x = np.linspace(min(xx), max(xx), 100)
    spline = UnivariateSpline(x, curve-np.max(curve)/2, s=0) # Find FWHM location
    r1, r2 = spline.roots()


def plot_violin(mpgs,
                mstar_cut_hard = 5e9,
                massive_cut=1e10,
                fname="",
                wdir='./',
                use_seaborn=True,
                scale = "width",
                pallette="muted",
                linewidth = 0.8,
                bw = 0.1,
                gridsize=10):
    
    dlt_all=[]
    dlo_all=[]
    dlM_all=[]
    dlm_all=[]
    mass_all=[]
    for igal, gal in enumerate(mpgs):
        if gal.data["mstar"][0] > mstar_cut_hard:
            dlt_all.append(gal.dlt) # Total
            dlo_all.append(gal.dlo) # Other
            dlM_all.append(gal.dlM) # Major
            dlm_all.append(gal.dlm) # minor
            mass_all.append(gal.data["mstar"][0]) # 
            #print(gal.dlM, gal.dlm, gal.dlt, gal.data["mstar"][0])
#        try:
#            mass_all.append(gal.data["mstar"][0].repeat(len(gal.dlm)))
#        except:

    
    dlM_all = remove_nan(np.array(dlM_all))
    dlm_all = remove_nan(np.array(dlm_all))
    dlo_all = remove_nan(np.array(dlo_all))
    mass_all = np.tile(remove_nan(np.array(mass_all)),3)
    #print(dlM_all)
    #print(dlm_all)
    #print(dlt_all)
    #print(mass_all)
    #x = x[~numpy.isnan(x)]

    if use_seaborn:
        import pandas as pd
        # Prepare data in Pandas DataFrame format. 
        val = np.concatenate((dlM_all, dlm_all, dlt_all))
        merger = np.concatenate((np.full(len(dlM_all), 0),
                                 np.full(len(dlm_all), 1),
                                 np.full(len(dlo_all), 2)))

        i_detected = val!=0
        massive = mass_all > massive_cut

        all_data = pd.DataFrame({"delta":val,
                                 "merger":merger,
                                 "massive":massive})
        #print(all_data["delta"], mass_all[massive])

        detected = pd.DataFrame({"delta":val[i_detected]
                                ,"merger":merger[[i_detected]]
                                ,"massive":mass_all[i_detected] > massive_cut})
        #print(sum((detected["merger"] == 0) * (detected["massive"])))
        
        fig, axs = plt.subplots(2,2,
                            sharex=True,
                            sharey=True)
        fig.set_size_inches(8,6)
        axs = axs.ravel()
        
        pos = [0,1,2]
        import seaborn.apionly as sns
        sns.set_style("whitegrid")
        #sns.set(font_scale=1.0)
        
        if violin:
            # All
            sns.violinplot(x="merger", y="delta",
                           data=all_data, palette=pallette, bw =bw,
                           scale=scale, ax=axs[0], inner="quartile",
                           linewidth=linewidth, gridsize=gridsize)
            axs[0].text(0, -0.7, str(sum(all_data["merger"] == 0))
                       , horizontalalignment="center")
            axs[0].text(1, -0.7, str(sum(all_data["merger"] == 1))
                        , horizontalalignment="center")
            axs[0].text(2, -0.7, str(sum(all_data["merger"] == 2))
                        , horizontalalignment="center")
            # Mass
            sns.violinplot(x="merger", y="delta",
                           data=all_data, hue="massive", palette=pallette,
                           bw =bw, scale=scale, split=True, ax=axs[1],
                           inner="quartile", linewidth=linewidth,
                           gridsize=gridsize)
            axs[1].text(0, -0.7, str(sum((all_data["merger"] == 0) * (~all_data["massive"]))) + \
                         " / " + str(sum((all_data["merger"] == 0) * (all_data["massive"])))
                       , horizontalalignment="center")
            axs[1].text(1, -0.7, str(sum((all_data["merger"] == 1) * (~all_data["massive"]))) + \
                         " / " + str(sum((all_data["merger"] == 1) * (all_data["massive"])))
                       , horizontalalignment="center")
            axs[1].text(2, -0.7, str(sum((all_data["merger"] == 2) * (~all_data["massive"]))) + \
                         " / " + str(sum((all_data["merger"] == 2) * (all_data["massive"])))
                       , horizontalalignment="center")
            # detected - all
            sns.violinplot(x="merger", y="delta",
                           data=detected, palette=pallette, bw =bw,
                           scale=scale, ax=axs[2], inner="quartile",
                           linewidth=linewidth, gridsize=gridsize)
            axs[2].text(0, -0.7, str(sum(detected["merger"] == 0))
                        , horizontalalignment="center")
            axs[2].text(1, -0.7, str(sum(detected["merger"] == 1))
                        , horizontalalignment="center")
            axs[2].text(2, -0.7, str(sum(detected["merger"] == 2))
                        , horizontalalignment="center")
            
            # detected - Mass
            last_ax = sns.violinplot(x="merger", y="delta",
                           data=detected, hue="massive", palette=pallette,
                           bw =bw, scale=scale, split=True, ax=axs[3],
                           inner="quartile", linewidth=linewidth, gridsize=gridsize)
            axs[3].text(0, -0.7, str(sum((detected["merger"] == 0) * (~detected["massive"]))) + \
                         " / " + str(sum((detected["merger"] == 0) * (detected["massive"])))
                        , horizontalalignment="center")
            axs[3].text(1, -0.7, str(sum((detected["merger"] == 1) * (~detected["massive"]))) + \
                         " / " + str(sum((detected["merger"] == 1) * (detected["massive"])))
                        , horizontalalignment="center")
            axs[3].text(2, -0.7, str(sum((detected["merger"] == 2) * (~detected["massive"]))) + \
                         " / " + str(sum((detected["merger"] == 2) * (detected["massive"])))
                        , horizontalalignment="center")
            # Remove Seaborn labels
            for ax in axs:
                ax.set_xlabel("")
                ax.set_ylabel("")
                ax.set_ylim([-0.8, 0.8])

            ### add custom labels

            # second row
            for ax in axs[2:4]:
                ax.set_xticks(pos)
                ax.set_xticklabels(labels = ["Major", "Minor", "non-merger"])
                ax.tick_params(labelsize=14)

            # left column
            for ax in [axs[0], axs[2]]:
                ax.tick_params(labelsize=14)
                ax.set_ylabel(r"$ \Delta \lambda_{R_{eff}}$", fontsize=16)   

            # right column
            for ax in [axs[1], axs[3]]:
            # Set legend #
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles, [r"$log_{10}M_{*} < $ " +"{:.1f}".format(np.log10(massive_cut)),
                                    r"$log_{10}M_{*} > $ " +"{:.1f}".format(np.log10(massive_cut))],
                                    loc='upper left',
                                    fontsize=11)
        else:
            ax = sns.swarmplot(data = data,
                               size=1)
        #ax.annotate(r"$ \log_{10}(M_{*}) >$" + "{:.1f}".format(\
        #                    np.log10(mstar_limit)), xy=[0.0,0.4],
        #                   fontsize=22)
    else:
        pos = [1,2,3]
        plt.violinplot([dlM_all, dlm_all, dlt_all], pos, 
                      points=30,
                      widths=0.5,
                      showmeans=True,
                      showextrema=True,
                      showmedians=True,
                      bw_method="silverman")
        
        #ax.annotate(r"$ \log_{10}(M_{*}) >$" + "{:.1f}".format(np.log10(mstar_limit)), xy=[0.7,0.5])
    plt.tight_layout()
    plt.savefig(fname, dpi=200)
    plt.close()
    print("Total {} galaxies above the mass cut {}".format(len(dlM_all), 
                                                           np.log10(mstar_cut_hard)))
    #plt.show()
    return last_ax


def plot_violin2(mpgs,
                mstar_cut_hard = 5e9,
                massive_cut=1e10,
                fname="",
                wdir='./',
                use_seaborn=True,
                scale = "width",
                pallette="muted",
                linewidth = 0.8,
                bw = 0.1,
                gridsize=10,
                violin = True):
    
    from matplotlib.ticker import NullFormatter
    
    # Per event data
    dl = []
    dm = []
    mr = []
    mass = []
    for gal in mpgs:
        if gal.data["mstar"][0] > mstar_cut_hard:
            if hasattr(gal, "merger"):
                if gal.merger is not None:
                    dl.extend(gal.merger.delta_l)
                    dm.extend(gal.merger.delta_m)
                    mr.extend(gal.merger.mr)
                    mass.extend([gal.data["mstar"][0]]*len(gal.merger.delta_l))
    
    i_ok = (np.array(mr) < 20 ) * (np.array(dl) > -1)
    dl = np.array(dl)[i_ok]
    dm = np.array(dm)[i_ok]
    mr = np.array(mr)[i_ok]
    mass = np.array(mass)[i_ok]
        

    # Prepare data in Pandas DataFrame format. 
    import pandas as pd
    
    # Per galaxy data
    dlt_all=[]
    dlo_all=[]
    dlM_all=[]
    dlm_all=[]
    mass_all=[]
    for igal, gal in enumerate(mpgs):
        if gal.data["mstar"][0] > mstar_cut_hard:
            dlt_all.append(gal.dlt) # Total
            dlo_all.append(gal.dlo) # Other
            dlM_all.append(gal.dlM) # Major
            dlm_all.append(gal.dlm) # minor
            mass_all.append(gal.data["mstar"][0]) # 

    dlM_all = remove_nan(np.array(dlM_all))
    dlm_all = remove_nan(np.array(dlm_all))
    dlo_all = remove_nan(np.array(dlo_all))
    mass_all = remove_nan(np.array(mass_all))
    #mass_all = np.tile(remove_nan(np.array(mass_all)),3)

    mr_major = 4
    dl_major = dl[mr < mr_major]
    dl_minor = dl[mr > mr_major]

    yll, ylr = (-0.5, 0.4)
    
    if use_seaborn:
        pallette="Paired"
        #sns.set_style("whitegrid", {'axes.grid' : False})
        plt.figure(1, figsize=(6,4))
        axHist = plt.axes([0,0,0.3,1])
        axViolin = plt.axes([0.31,0,1,1])
        
        sns.distplot(dl_major, color="r",hist=True
                     , kde_kws={"shade": True}, vertical=True
                     , ax=axHist, label="major mergers")
        
        sns.distplot(dl_minor, color="b",hist=True
                     , kde_kws={"shade": True}, vertical=True
                     , ax=axHist, label="minor mergers")
        
        sns.distplot(dlo_all, color="g",hist=True
                     , kde_kws={"shade": True}, vertical=True
                     , ax=axHist, label="non-mergers")
        
        axHist.xaxis.set_major_formatter(NullFormatter())
        axHist.set_ylabel(r"$\Delta \lambda_{R_{eff}}$")
        axHist.set_xlabel("Normalized probability")
        axHist.set_ylim([-0.7,0.6])
        
        axHist.legend(fontsize=12)
        axHist.yaxis.grid()
        
        ## Violin plot data
        Merger_alone=False
        if Merger_alone:
            all_data = pd.DataFrame({"delta":dl,
                                 "merger":mr < mr_major,
                                 "massive":mass > massive_cut})
        else:
            
            mass_concat = np.concatenate((mass, mass_all))
            mr_tag = np.concatenate((np.full(len(mr), 1), np.full(len(dlo_all), 2)))
            mr_tag[mr < mr_major] = 0
            is_massive = np.concatenate((mass, mass_all)) > massive_cut
            all_data = pd.DataFrame({"delta":np.concatenate((dl,dlo_all)),
                                 "merger":mr_tag,
                                 "massive":is_massive})
        
            sns.violinplot(x="merger", y="delta",
                           data=all_data, hue="massive", palette=pallette,
                           bw=bw, scale=scale, split=True, ax=axViolin,
                           inner="quartile",
                           linewidth=linewidth,
                           gridsize=gridsize,
                           )
        axViolin.text(0, -0.65, str(sum((mr_tag == 0) * ~is_massive)) + "/" +\
                      str(sum((mr_tag == 0) * is_massive)), horizontalalignment="center")
        axViolin.text(1, -0.65, str(sum((mr_tag == 1) * ~is_massive)) + "/" +\
                      str(sum((mr_tag == 1) * is_massive))
                      , horizontalalignment="center")
        axViolin.text(2, -0.65, str(sum((mr_tag == 2) * ~is_massive)) + "/" +\
                      str(sum((mr_tag == 2) * is_massive))
                      , horizontalalignment="center")    
        
        axViolin.set_xticks([0,1,2])
        axViolin.set_xticklabels(labels = ["Major", "Minor", "non-merger"])
        axViolin.tick_params(labelsize=14)
        axViolin.yaxis.set_major_formatter(NullFormatter())
        axViolin.set_ylabel("")
        axViolin.set_xlabel("")
        axViolin.set_ylim([-0.7,0.6])
        
        handles, labels = axViolin.get_legend_handles_labels()
        axViolin.legend(handles, [r"$log_{10}M_{*} < $ " +"{:.1f}".format(np.log10(massive_cut)),
                                  r"$log_{10}M_{*} > $ " +"{:.1f}".format(np.log10(massive_cut))],
                            loc='upper left',
                            fontsize=12)
        axViolin.yaxis.grid()
        
    else:
        same_x = True
        # same x scale 
        if same_x:
            hist_major, bin_edges = np.histogram(dl_major, bins=10, range=(yll,ylr))
            hist_minor, bin_edges = np.histogram(dl_minor, bins=10, range=(yll,ylr))

            mM = max(hist_major)
            mm = max(hist_minor)

            width_M = mM/(mM + mm) * 0.4
            width_m = mm/(mM + mm) * 0.4
        else:
            width_M = 0.2
            width_m = 0.2

        # definitions for the axes
        left, width = 0.1, 0.4
        bottom, height = 0.1, 0.7

        #rect_scatter = [left, bottom, width, height]
        #rect_major = [left + width, bottom, width_M, height]
        #rect_minor = [left + width + width_M, bottom, width_m, height]

        # start with a rectangular Figure
        plt.figure(1, figsize=(16,8))

        #axScatter = plt.axes(rect_scatter)
        #axHist1 = plt.axes(rect_major)
        #axHist2 = plt.axes(rect_minor)
        fig, axHist2 = plt.subplots()

        # no labels
        nullfmt   = NullFormatter() # no labels
        #axHist1.yaxis.set_major_formatter(nullfmt)
        axHist2.yaxis.set_major_formatter(nullfmt)


        hm = axHist2.hist(dl_minor, orientation='horizontal', color='green', range=(yll,ylr))
        gaussian_fit(axHist2, dl_minor, hm[1][1] - hm[1][0], color='greenyellow')
        axHist2.set_xlim((0, 1.1*max(hm[0])))

        hM = axHist2.hist(dl_major, orientation='horizontal', color='blue', range=(yll,ylr))
        gaussian_fit(axHist2, dl_major, hM[1][1] - hM[1][0], color='cyan')


        #axHist2.barh(np.arange(len(hist_minor)), hist_minor, color='green')


        #axHist1.set_ylim(axScatter.get_ylim())
        axHist2.set_xlabel("# of mergers", x=0)
    #fig.tight_layout()
    #plt.show()
    plt.savefig(fname + ".png", dpi=200, bbox_inches="tight")
    plt.savefig(fname + ".pdf", bbox_inches='tight') # eps does NOT support transparency!
    plt.close()
    return av


def plot_simple(mpgs, wdir='./', suffix=""):
    fig, ax = plt.subplots(2)
    for igal, gal in enumerate(mpgs):
        if gal.merger is not None:
            for mr, xx, delta in zip(gal.merger.mr, gal.merger.nout, gal.merger.delta_l):
                ax[0].scatter(mr, delta)
        ax[1].scatter([1],gal.dlM) 
        ax[1].scatter([3],gal.dlm) 
        ax[1].scatter([5],gal.dlo)
    #plt.show()
    ax[0].set_xscale("log")
    plt.savefig(wdir + "mma_simple" + suffix + ".png")
    plt.close()


# In[3]:

def plot_hist(mpgs, wdir='./', suffix=""):
    all_mr = []
    all_delta =[]
    fig, ax = plt.subplots()
    for gal in mpgs:
        if gal.merger is not None:
            for mr, delta in zip(gal.merger.mr, gal.merger.delta_l):
                all_mr.append(mr)
                all_delta.append(delta)

    all_mr = np.array(all_mr)
    all_delta = np.array(all_delta)
    plt.ioff()
    ax.hist(all_delta[all_mr < 3], alpha=0.3, range=[-0.7, 0.5],
             facecolor='none', edgecolor='r', linestyle='dashed',
             lw=3)
    ax.hist(all_delta[(all_mr < 10) * (all_mr > 3)], alpha=0.3, range=[-0.7, 0.5],
             facecolor='none', edgecolor='g', linestyle='solid', lw=3)
    ax.hist(all_delta[all_mr > 10], alpha=0.3, range=[-0.7, 0.5],
             facecolor='none', edgecolor='b', linestyle="-.", lw=3)
    plt.savefig(wdir + 'MMA_hist' + suffix +'.png')
    plt.close()
    
    print(len(all_mr))



# In[4]:

def filter_small_mergers(mm, window=7):
    """
    If multiple mergers occur within short period of time, leave the largest.
    """
    if mm is None:
        return
    
    if len(mm.nout) == 1:
        return

    neighbor = np.zeros(len(mm.mr))
    # Tag mergers with neighbor
    delimeter = 1
    for i in range(len(mm.mr)-1):
        #print(mm.nout[i], mm.nout[i+1], delimeter)
        neighbor[i]=delimeter
        if abs(mm.nout[i] - mm.nout[i+1]) > window:
            delimeter +=1

        
        #print(neighbor)
    # last element
    neighbor[-1] = delimeter
    #if abs(mm.nout[-2] - mm.nout[-1]) < window:
    #    neighbor[-1] = delimeter
    #else:
    #    neighbor[-1] = delimeter + 1

    #print(mm.nout)
    #print(neighbor)

    ind_ok = []
    # neighbor may have 0,1, or may start from 2.
    for imerger in range(delimeter+1):
        ii = np.where(neighbor == imerger)[0]
        if len(ii) > 0:
            ind_ok.append(ii[np.argmin(mm.mr[ii])])

    #print("finally", neighbor[ind_ok],mm.mr, mm.mr[ind_ok])
    mm.nout = mm.nout[ind_ok]
    mm.mr = mm.mr[ind_ok]
    mm.nout_ini = mm.nout_ini[ind_ok]

    
def remove_nan(x):
    return x[~np.isnan(x)]


def l_at_smoothed_r(gal, npix_per_reff=5):
    """
    Lambda_r measured at 4 or 5?
    I used 5 in calculating lambda. It's an average over from 1.0Reff to 1.2Reff, may be?
    """
    n_valid_points = sum(gal.data["reff"] > 0)
    new_l_arr = np.zeros(n_valid_points)
    new_reff = mma.smooth(gal.data["reff"])#[gal.data["reff"] > 0])
    for i in range(n_valid_points):
        try:
            lambdar = gal.data["lambda_arr"][i]
            ind_org = npix_per_reff
            #i_new =  # 0-indexed.
            ind_new = new_reff[i]/gal.data["reff"][i] * npix_per_reff
            il = np.fix(ind_new).astype(int)
            ir = il + 1
            if ir >= len(lambdar):
                new_l_arr[i] = lambdar[-1]
            else:
                new_l_arr[i] = lambdar[il]*(ir-ind_new) + lambdar[ir]*(ind_new-il)
                #print(lambdar[ind_org] / new_l_arr[i])
                #new_l_arr[i] = lambdar[il]*(ir-ind_org) + lambdar[ir]*(ind_org-il)
        except:
            new_l_arr[i] = gal.data["lambda_arr"][i][ind_org] # = 0 with bad measurements.
            
    return new_l_arr


def individual_lambda_evol(mpgs, fname="figs/individual_lambda_evol.pdf"):
    from matplotlib.backends.backend_pdf import PdfPages
    fig, ax = plt.subplots(1, sharex=True)
    with PdfPages(fname) as pdf:
        for gal in mpgs:
            ax.plot(gal.smoothed_lambda)
            ax.set_title(str(gal.cluster) + "  " + str(gal.idxs[0]))
            pdf.savefig()
            ax.clear()
            
            

            
def find_merger_epochs(alltrees, 
                       idx_all, 
                       mpgs, 
                       nout_ini=37, 
                       dist_gal_scale_in=2.0,
                       dist_gal_scale_out =3.0,
                       min_mass_ratio = 0.05,
                       mass_ratio="early",
                       verbose=False,
                       do_plot = False,
                       pdf_fname='./merger_ratio.pdf',
                       max_rgal=40):
    """
    Parameters
    ----------
    dist_gal_scale 
        if two galaxies are closer than dist_gal_scale * (sum of raidus of the two),
        that epoch is the nout_init_merger.
    nout_ini
        blabla
    """
    gal_list=[]
    mr_list=[]
    nout_list=[]
    nout_ini_list=[] # initial time when two halos(Galaxy stellar components in this case) overlap. 

    #for idx in idx_all:
    if do_plot:
        from matplotlib.backends.backend_pdf import PdfPages
        fig, ax = plt.subplots(1, sharex=True)
        pdf = PdfPages(pdf_fname)
    else:
        ax = None

    for gal in mpgs:
        idx = gal.data['idx'][0]
        # full tree of a galaxy
        atree = ctu.extract_a_tree(alltrees.data, idx)

        # main progenitor tree
        main = ctu.extract_main_tree(atree, idx)

        main_nout = main['nout'].flatten()
        i_nout_ok = main_nout > nout_ini
        main = main[i_nout_ok]
        #x_nout = x_nout[i_nout_ok]
        pos = np.zeros((3,len(main)))
        pos[0,:] = main['x']
        pos[1,:] = main['y']
        pos[2,:] = main['z']

        
        ## 왜 len(main)?
        mass_ratios_single = np.zeros(len(main))
        nout_inits = np.zeros(len(main))
        #print("log M* ={}".format(np.log10(gal.data["mstar"][0])))
        
        ## Substitute main["r"] with gal.data["reff"]
        if len(main) < len(gal.smoothed_r):
            main["r"] = gal.smoothed_r[:len(main)]
        elif len(main) > len(gal.smoothed_r):
            main["r"][:len(gal.smoothed_r)] = gal.smoothed_r
        elif len(main) == len(gal.smoothed_r):
            main["r"] = gal.smoothed_r
        
        for i, nout in enumerate(main['nout']):
            # merger ratio
            # First, check if there are multiple progenitors.
            i_prgs = np.where(atree['desc_id'] == main['id'][i])[0]
            
            # multiple prgs = merger
            if len(i_prgs) > 1:
                if verbose: 
                    print("idx:{}, {} Progenitors at nout = {}".format(main["id"][0], len(i_prgs), nout))
                
                #print("i, inout_mpgs", i, np.where(gal.nouts == nout)[0])
                #print(gal.data["reff"][i])
                
                # Mass ratio must be calculated inside get_merger_info. 
                
                id_prgs = atree['id'][i_prgs]
                mass_prgs = atree['m'][i_prgs]
                #m_r = mass_prgs / max(mass_prgs)
                
                # Progenitor with maximum mass at the FINAL COALESCENCE is the main progenitor.
                # Others are satellites.
                sats = id_prgs[mass_prgs < max(mass_prgs)]
                
                mass_ratios_now=[]
                nout_inits_now=[]
                
                # loop over satellites at a given nout.
                for this_sat in sats:
                    n_i_t, mass_this_sat = get_merger_info(main, atree, this_sat,
                                                           dist_gal_scale_in=dist_gal_scale_in,
                                                           dist_gal_scale_out = dist_gal_scale_out,
                                                           do_plot=do_plot,
                                                           ax=ax, 
                                                           max_rgal=max_rgal)
                    
                    mass_ratio = mass_this_sat / max(mass_prgs)
                    if do_plot:
                        ax.text(40, 600, "{:.3f}".format(mass_ratio))
                        pdf.savefig()
                        ax.clear()
                    if mass_ratio > min_mass_ratio:
                        nout_inits_now.append(n_i_t)
                        mass_ratios_now.append(mass_ratio)
                if len(nout_inits_now) > 0:
                    i_main_merger = np.argmax(np.array(mass_ratios_now))
                    mass_ratios_single[i] = mass_ratios_now[i_main_merger]
                    nout_inits[i] = nout_inits_now[i_main_merger]
                else:
                    mass_ratios_single[i] = -1
                    nout_inits[i] = -1

                #if verbose:
                #    print(" Mass ratios {} at nout = {}: ".format(m_r, nout_inits[i]))
                #if do_plot:
                #    pdf.savefig()
                    #ax.clear()

            else:
                mass_ratios_single[i] = 0
            ##--------------------------------------------------
                    

        ind_ok = np.where(mass_ratios_single > min_mass_ratio)[0]
        if len(ind_ok) > 0:
            # if a satellite oscillates around the host, 
            # it could be identified as multiple mergers with short time interval. 
            # leave only the first passage / merger.
            # No, it doesn't happen in ConsistentTrees.

            #good =[]
            #for i in range(len(ind_ok)-1):
            #    if ind_ok[i+1] > ind_ok[i] + 2:
            #        good.append(ind_ok[i])
            #good.append(ind_ok[-1])
            #ind_ok = good
            mr = 1./mass_ratios_single[ind_ok]

            gal_list.append(idx)
            mr_list.append(mr)
            nout_list.append(main_nout[ind_ok])    
            nout_ini_list.append(nout_inits[ind_ok])
            #print(idx)
    if do_plot:
        pdf.close()

    inds=[]
    for i, gal in enumerate(mpgs):
        galid = gal.data['idx'][0]
        ind = np.where(galid == gal_list)[0]
        if len(ind) > 0:
            inds.append(i)
            merger = Merger()
            merger.mr = mr_list[ind]
            merger.nout = nout_list[ind]
            merger.nout_ini = nout_ini_list[ind]
            gal.merger = merger
        else:
            gal.merger = None
    
    
class Merger():
    pass


def get_merger_info(main, atree, sat_root_idx,
                    dist_gal_scale_in=1.0,
                    dist_gal_scale_out=2.0,
                    do_plot=True,
                    ax=None,
                    max_rgal=40):
    """
    Returns merger mass ratio at the BEGINNING of the merger. 
    nout_init_this_merger, mass_this_merger = get_merger_info()
    
    Main["r"] = gal.smoothed_r.
    So I can use the main galaxy radius to decide whether a satellite is
    touching/affecting the main galaxy or not.
    """
    satellite = ctu.extract_main_tree(atree, sat_root_idx, no_subset=True)
    
    # nouts where both of galaxies' trees exist.
    nout_min = max([min(main['nout']), min(satellite['nout'])])
    idx_main = main["id"][0]
    i_main_ok = (main['nout'] >= nout_min) * (main['nout'] <= max(satellite["nout"]))
    i_sat_ok = (satellite['nout'] >= nout_min)
    satellite = satellite[i_sat_ok]
    main_this=main[i_main_ok]

    # distances at all valid nouts in kpc.
    dd = np.sqrt(np.square(main_this["x"] - satellite['x']) \
               + np.square(main_this["y"] - satellite['y']) \
               + np.square(main_this["z"] - satellite['z'])) * 1e3
    #print(dd)
    rgal_tot = main_this['r']# + satellite['rvir'])
    rgal_tot[rgal_tot > max_rgal] = max_rgal
    #print(rgal_tot)
    #print(main_this["nout"])
    #print(" Galaxy sizes : main {}, and the second {}, and the sum {}".format(
    #        main['r'][i_main_ok], satellite['r'], rgal_tot))
    #print(" dd :", dd)
    if do_plot:
        ax.plot(main_this['nout'], rgal_tot, label="main")
        ax = plt.gca()
        ax.plot(main_this['nout'], rgal_tot * dist_gal_scale_in, label="threshold in")
        ax.plot(main_this['nout'], rgal_tot * dist_gal_scale_out, label="threshold out")
        ax.plot(satellite['nout'], dd, label="distance")
        
        #plt.set_yscale('log')
        
        ax.set_ylim([1,1e3])
        ax.set_xlim([30,200])
        ax.set_title("main:" + str(idx_main) +  " sat:" + str(sat_root_idx))
        ax.legend()

    
    # A steep plunging orbit may cross the threshold in one snapshot.
    # in such cases, the final snapshot should be within the 'out' threshold,
    # and then the beginning of merger is thought to be the one right before it merge.
    # So any galaxy that entered closer than the dis_gal_scale_OUT are taken into account.
    # This assumes that these galaxies DO merge at the next snapshot. 
    # So make sure that the tree is correct.
    if sum(dist_gal_scale_out * rgal_tot > dd) > 0:
        # First close encounter is technically the beginning of merger,
        # but in practice that could be merely a flyby, 
        # and whether they will merger soon or not is not known. 
        # I can't call an encounter a merger if the encounter will end up merging in 100Gyrs.        
        #nout_init_this = min(satellite['nout'][dist_gal_scale * rgal_tot < dd])

        i_dist_out = np.where(dist_gal_scale_out * rgal_tot < dd)[0]
        i_dist_bet = np.where((dist_gal_scale_out * rgal_tot > dd) \
                              * dist_gal_scale_in * rgal_tot < dd)[0]
        i_dist_in = np.where(dist_gal_scale_in * rgal_tot >= dd)[0]
        if len(i_dist_in) > 0:
            if len(i_dist_out) > 0:
                # index 0 = last snapshot.
                i_dist_out_last = min(i_dist_out)
                # earliest time since when the galaxy doesn't go further than the threshold.
                i_good = i_dist_in < i_dist_out_last
                if sum(i_good) > 1:
                    i_dist_final = max(i_dist_in[i_good]) + 1 # max = earliest.
                    nout_init_this = satellite['nout'][i_dist_final]
                    mass_this = satellite['m'][satellite['nout'] == nout_init_this].squeeze()
                    if do_plot: ax.text(40, 800, "good 1")
                elif sum(i_good) == 1:
                    i_dist_final = i_dist_in[i_good] + 1
                    if do_plot: ax.text(40, 800, "good 2")
                    nout_init_this = satellite['nout'][i_dist_final]
                    mass_this = satellite['m'][satellite['nout'] == nout_init_this].squeeze()
                else:
                    #print(i_dist_in, i_dist_out_last, "What's going on???")
                    if do_plot: ax.text(40, 800, "unphysical...")
                    i_dist_final = 0 # dummy
                    nout_init_this = -1
                    mass_this = 0
            else:
                if do_plot: ax.text(40, 800, "Warning.. No i_dist_out")
                i_dist_final = max(i_dist_in)
                nout_init_this = satellite['nout'][i_dist_final]
                mass_this = satellite['m'][satellite['nout'] == nout_init_this].squeeze()
            #if len(i_dist_final) > 0:
            
        elif 0 in i_dist_bet:
            if do_plot: ax.text(40, 800, "bet")
            # A steep plunging orbit may cross the threshold in one snapshot.
            # in such cases, the final snapshot should be within the 'out' threshold,
            # and then the beginning of merger is thought to be the one right before it merge.
            i_dist_final = 0
            nout_init_this = satellite['nout'][0]
            mass_this = satellite['m'][0].squeeze()
        else:
            if do_plot: ax.text(40, 800, "bad 1")
            i_dist_final = -1
            nout_init_this = -1
        #print(nout_init_this, dd[i_dist_final])
        if do_plot: ax.scatter(nout_init_this, dd[i_dist_final])
        
    else:
        #print("None")
        nout_init_this = -1
        mass_this = 0
    
    return nout_init_this, mass_this




def measure_delta(mpgs, 
                  nout_ini=37, 
                  nout_fi = 187,
                  wdir='./',
                  dt_after = 0.5,
                  dt_before = 0.5,
                  dt_settle = 0.5,
                  savefig=False,
                  figname="figs/measure_delta"):
    """
    Measure lambda change at every merger and draw lambda evolution with merger events.
    time span over which average lambda value is measured is given in Gyr unit.
    
    
    lambda_before = average(lambda[nout_begin-dt_before] ~ lambda[nout_begin])
    lambda_after  = average(lambda[nout_merger-dt_before] ~ lambda[nout_merger])
    ** nout_merger = nout at final coalescence.
    
    Parameters
    ----------
    dt_after:
        d
    
    dt_before:
        ddd
    
    dt_settle:
        aa


    galaxies must have smoothed quantaties.
    """
    
    from utils.util import dgyr2dnout    
    if savefig:
        from matplotlib.backends.backend_pdf import PdfPages
        fig, ax = plt.subplots(1,2, sharex=True)
        plt.subplots_adjust(hspace=0.001)
        fig.set_size_inches(8,4)    
        pdf = PdfPages(wdir + figname + '.pdf')
        
    for gal in mpgs:
        ind_nout = gal.nouts > nout_ini
        gal.nouts = gal.nouts[ind_nout]
        gal.data = gal.data[ind_nout]
        if savefig:
            ax[0].plot(gal.nouts, np.log10(gal.data['mstar']), 'b:')
            ax[0].set_xlim([50,190])
            plt.suptitle(str(gal.cluster) + " - " + 
                         str(gal.ids[0]) + ", " + 
                         str(gal.idxs[0]) + 
                         "  {:.2e}".format(gal.data["mstar"][0]))
            ax[1].plot(gal.nouts, gal.smoothed_lambda, 'black')

        if gal.merger is not None:
            delta_lambda =[]
            delta_mass = []
            n_after = []
            n_before = []
            
            for mr, nout_merger, nout_begin in zip(gal.merger.mr, gal.merger.nout, gal.merger.nout_ini):
                # nout_ini가 너무 옛날인 경우 gal.nouts에 nout_ini - dt가 없음. 
                nout_merger =min([nout_merger + 3, 187])
                if True:
                    # index of merger (final coalescence)
                    if nout_merger < min(gal.nouts):
                        print("Too early merger")
                        # Too early merger
                        delta_lambda.append(-1)
                        delta_mass.append(-1)
                        n_after.append(0)
                        n_before.append(0)
                        continue
                    
                    i_nout = np.where(gal.nouts == nout_merger)[0]                    
                    # index of the begining of merger (determined in find_merger_epoch())
                    iini_nout = np.where(gal.nouts == nout_begin)[0]
                    if len(iini_nout) == 0:
                        # Tree may be longer than valid lambda measurement points.
                        print("end of merger in gal.nouts, but begining of merger is not.")
                        iini_nout = len(gal.nouts) -1

                    # snapshot number apart by dt from nout
                    nout_settle = dgyr2dnout(dt_settle, nout_merger)
                    nout_after = dgyr2dnout(dt_settle + dt_after, nout_merger)
                    
                    # indices of the snapshots.
                    #if min(gal.nouts) < nout_after:
                    if min(gal.nouts) < nout_settle:
                        inout_after = np.where(gal.nouts == nout_after)[0]
                    else:
                        delta_lambda.append(-1)
                        delta_mass.append(-1)
                        n_after.append(0)
                        n_before.append(0)
                        print("end of merger + dt_settle is earlier than the min nout.")
                        continue                       
                        
                    if i_nout == 0:
                        l_inds_after = [0]
                    else:
                        inout_settle = np.where(gal.nouts == nout_settle)[0]
                        #print("inout_after, i_nout, xx, gal.nouts", inout_after, i_nout, xx, gal.nouts)
                        #l_inds_after = range(max(inout_after), i_nout)
                        l_inds_after = range(inout_after, min([i_nout, inout_settle])) # index 0 = nout 187.
                    
                    # quantities at measurement points
                    nouts_after = gal.nouts[l_inds_after]
                    l_after = gal.smoothed_lambda[l_inds_after]
                    m_after = gal.data['mstar'][l_inds_after]
                    lambda_after = np.average(l_after)
                    mass_after = np.average(m_after)
                        
                        
                    nout_before = dgyr2dnout(-1 * dt_before, nout_begin)
                    nout_before = max([nout_before, min(gal.nouts)])
                    # Tree may be longer than valid lambda measurement points.
                    inout_before = np.where(gal.nouts == nout_before)[0]

                    # if time_{merger} + dt > t_{z=0}
                    #l_inds_before = range(iini_nout,min([np.where(gal.data["reff"] > 0)[0][-1], inout_before]))
                    #print("x2, iini_nout, inout_before", x2, iini_nout, inout_before)
                    l_inds_before = range(iini_nout, inout_before + 1)
                    nouts_before = gal.nouts[l_inds_before]
                    l_before = gal.smoothed_lambda[l_inds_before]
                    m_before = gal.data['mstar'][l_inds_before]
                    lambda_before = np.average(l_before)
                    mass_before = np.average(m_before)

                    delta_lambda.append(lambda_after - lambda_before)
                    delta_mass.append(mass_after - mass_before)
                    n_after.append(nout_after)
                    n_before.append(nout_before)
                    if savefig:
                        ax[1].plot(nouts_after,l_after, 'g-')                                           
                        # Check again.
                        #nn = range(
                        # Average value
                        ax[1].plot(nouts_after, [lambda_after]*len(nouts_after), "g:")                    

                        ax[1].plot(nouts_before,l_before, 'r-')
                        #nn = range(max([nout_ini, min(nouts_before) - 5]), max(nouts_before) + 5)
                        ax[1].plot(nouts_before, [lambda_before]*len(nouts_before), "r:")

                        ax[0].axvline(nout_merger, linestyle=':')
                        ax[0].annotate("{:.1f}".format(mr), xy=(nout_merger,0.8))
                        ax[1].axvline(nout_merger, linestyle=':')
                        ax[1].annotate("{:.1f}".format(mr), xy=(nout_merger,0.8))
                        ax[1].axvline(nout_merger, linestyle=':')
                        #ax[1].axvline(x2, linestyle=':', c='g')
                #except:
                else:
                    delta_lambda.append(-1)
                    delta_mass.append(-1)
                    n_after.append(0)
                    n_before.append(0)


            gal.merger.delta_l = np.array(delta_lambda)
            gal.merger.delta_m = np.array(delta_mass)    
            gal.merger.nout_after = np.array(n_after)
            gal.merger.nout_before = np.array(n_before)
        if savefig:
            pdf.savefig()
            ax[0].clear()
            ax[1].clear()


    if savefig:
        pdf.close()




def Maj_min_acc_ratio(mpgs, dt=5, major_ratio=4):
    for igal, gal in enumerate(mpgs):
        delta_lambda_tot = np.average(gal.data['lambda_r'][:dt]) - np.average(gal.data['lambda_r'][-dt:])

        delta_lambda_major=0
        delta_lambda_minor=0
        if gal.merger is not None:

    #        for mr, xx, delta in zip(gal.merger.mr, gal.merger.nout, gal.merger.delta):
    #            ax[0].scatter(mr, delta)
            #print(gal.merger.mr, gal.merger.delta_l)
            i_major = np.where(gal.merger.mr <= major_ratio)[0]
            if len(i_major) > 0:
                #print(gal.merger.mr, gal.merger.delta_l, i_major, "M")
                delta_lambda_major = 0
                for iii in i_major:
                    if gal.merger.delta_l[iii] > -1:
                        delta_lambda_major += gal.merger.delta_l[iii]
                    #delta_lambda_major = np.sum(gal.merger.delta_l[i_major])

            i_minor = np.where(gal.merger.mr > major_ratio)[0]
            delta_lambda_minor = 0
            if len(i_minor) > 0:
                #print(gal.merger.delta_l, i_minor, "m")
                for iii in i_minor:
                    if gal.merger.delta_l[iii] > -1:
                        delta_lambda_minor += gal.merger.delta_l[iii]
#                    delta_lambda_minor = sum(gal.merger.delta_l[i_minor])


        delta_lambda_other = delta_lambda_tot - delta_lambda_major - delta_lambda_minor
        gal.dlt = delta_lambda_tot
        gal.dlo = delta_lambda_other
        gal.dlM = delta_lambda_major
        gal.dlm = delta_lambda_minor
