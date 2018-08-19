# coding: utf-8

# In[ ]:

# required : 
# info, galaxy catalog, halo catalog, tree, lambda catalog (catalog_187.pickle)
#
#
#
#
def mma_all(clusters, nout_fi=187, nout_ini=57, cdir='easy/',
            savefig_pdf=False, savefig_simple=True, savefig_violin=True):
    import utils.sampling as smp
    import matplotlib.pyplot as plt
    import tree
    import pickle
    import tree.halomodule as hmo
    import numpy as np
    import draw
    import load
    import analysis.evol_lambda as evl
    import analysis.Major_Minor_accretion as mma
    import tree.ctutils as ctu
    import utils.match as mtc

    #nout_fi = 187
    #nout_ini = 57
    #cdir = 'easy/'
    
    for wdir in clusters:   
        #wdir = './'
        dir_out = wdir 
        suffix = wdir[:-1]

        info = load.info.Info(base=wdir, nout=187, load=True)
        gg = hmo.Halo(nout=187, base=wdir, load=True, is_gal=True)
        hh = hmo.Halo(nout=187, base=wdir, load=True, is_gal=False)

        good_gals = mma.close_gals(hh, gg, rscale=3)

        alltrees = ctu.load_tree(wdir, is_gal=True)
        # Serialize catalogs. -> Only main galaxies
        # main galaxy list

        ad = alltrees.data
        tn = ad[ad['nout'] == nout_fi]

        cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout_fi) + '.pickle', 'rb'))
        # Only good galaxies are saved in the catalogue.
        # But there is no mass cut in good_gals list.
        # And match_list_ind doest not handle cases where the array a is shorter than array b.
        # Need another matching funciton that automatically excludes missing values.
        good_gals = np.intersect1d(cat['id'], good_gals)
        ind_ok = mtc.match_list_ind(cat['id'], good_gals, allow_swap=False)
        cat = cat[ind_ok]

        # Some galaxies are not in the tree. So can't follow their main progenitors.
        cat = cat[~np.isnan(cat['idx'])]
        idx_all = cat['idx'].astype(int) # fix so that idx is int from the beginning. 
        #print(good_gals)


        mpgs = mma.compile_mpgs(alltrees, idx_all, base=wdir, cdir=cdir, nout_ini=nout_ini, nout_fi=nout_fi)
        mma.find_merger_epochs(alltrees, idx_all, mpgs, nout_ini=nout_ini)
        mma.measure_delta_lambda(mpgs, dt=7, nout_ini=nout_ini)
        mma.Maj_min_acc_ratio(mpgs)
        pickle.dump(mpgs, open(dir_out + "mpgs" + suffix + ".pickle", "wb"))

        if savefig_pdf:
            from matplotlib.backends.backend_pdf import PdfPages

            fig, ax = plt.subplots(2,2, sharex=True)
            plt.subplots_adjust(hspace=0.001)

            with PdfPages(wdir + 'multipage_pdf.pdf') as pdf:
                for gal in mpgs:
                    ind_nout = gal.nouts > nout_ini
                    gal.nouts = gal.nouts[ind_nout]
                    gal.data = gal.data[ind_nout]
                    #smoothed_lambda = medfilt(gal.data['lambda_r'], kernel_size=5)

                    ax[0,0].scatter(gal.nouts, np.log10(gal.data['mstar']))
                    ax[0,0].set_xlim([50,190])
                    ax[0,0].set_title(str(gal.ids[0]) + ", " + str(gal.idxs[0]))
                    #ax[0].set_ylim([8,13])
                    ax[1,0].plot(gal.nouts, gal.smoothed_lambda, 'r-')
                    ylim_mid = ax[1,0].get_ylim() # match ylim of original and smoothed plot
                    ax[0,1].plot(gal.nouts, gal.data['lambda_r'], 'r-')
                    ax[0,1].set_ylim(ylim_mid)



                    if gal.merger is not None:
                        delta_lambda =[]
                        for mr, xx, x2 in zip(gal.merger.mr, gal.merger.nout, gal.merger.nout_ini):
                            ax[0,0].axvline(xx, linestyle=':')
                            ax[0,0].annotate("{:.1f}".format(mr), xy=(xx,0.8))
                            ax[0,1].axvline(xx, linestyle=':')
                            ax[0,1].annotate("{:.1f}".format(mr), xy=(xx,0.8))
                            ax[1,0].axvline(xx, linestyle=':')
                            ax[1,0].axvline(x2, linestyle=':', c='g')

                            i_nout = np.where(gal.nouts == xx)[0]
                            #print(max([0, i_nout - dt]))
                            #print(i_nout)

                            # nout[0] = 187
                            #lambda_after = np.average(smoothed_lambda[max([0, i_nout - dt]) : i_nout])
                            #lambda_before = np.average(smoothed_lambda[i_nout:min([len(gal.data), i_nout + dt])])
                            #delta_lambda.append(lambda_after - lambda_before)
                        #gal.merger.delta = np.array(delta_lambda)
                    pdf.savefig()
                    ax[0,0].clear()
                    ax[0,1].clear()
                    ax[1,0].clear()

            plt.close()

        if savefig_simple:
            fig, ax = plt.subplots(2)
            for igal, gal in enumerate(mpgs):
                if gal.merger is not None:
                    for mr, xx, delta in zip(gal.merger.mr, gal.merger.nout, gal.merger.delta):
                        ax[0].scatter(mr, delta)
                ax[1].scatter([1],gal.dlM) 
                ax[1].scatter([3],gal.dlm) 
                ax[1].scatter([5],gal.dlo)
            #plt.show()
            plt.savefig(dir_out + "mma_simple" + suffix + ".png")

        if savefig_violin:
            dlt_all=[]
            dlo_all=[]
            dlM_all=[]
            dlm_all=[]
            for igal, gal in enumerate(mpgs):
                dlt_all.append(gal.dlt)
                dlo_all.append(gal.dlo)
                dlM_all.append(gal.dlM)
                dlm_all.append(gal.dlm)


            data = [dlM_all, dlm_all, dlt_all]
            pos = [1,2,3]
            fig, ax = plt.subplots()
            ax.violinplot(data, pos, points=20, widths=0.3,
                                  showmeans=True, showextrema=True, showmedians=True)
            ax.set_ylim([-0.8, 0.3])
            ax.tick_params(labelsize=18)
            plt.savefig(dir_out + "mma_violin" + suffix + ".png")


# In[ ]:

def compile_and_plot(clusters, nout_fi=187, nout_ini=57, cdir='easy/', fout="mma_violin_all.png"):
    import pickle
    import matplotlib.pyplot as plt
    """
        compile mpgs.pickle of all clusters and draw a violing plot of all galaxies.
    """
    dlt_all=[]
    dlo_all=[]
    dlM_all=[]
    dlm_all=[]
    cluster_name=[]
    for wdir in clusters:
        dir_out = wdir 
        suffix = wdir[:-1]
        mpgs = pickle.load(open(dir_out + "mpgs" + suffix + ".pickle", "rb"))
        
        # compile
        for gal in mpgs:
            dlt_all.append(gal.dlt)
            dlo_all.append(gal.dlo)
            dlM_all.append(gal.dlM)
            dlm_all.append(gal.dlm)
            cluster_name.append(wdir)
            print(gal.dlm)
    print(len(dlM_all))
    data = [dlM_all, dlm_all, dlt_all]
    pos = [1,2,3]
    fig, ax = plt.subplots()
    ax.violinplot(data, pos, points=20, widths=0.3,
                          showmeans=True, showextrema=True, showmedians=True)
    ax.set_ylim([-0.8, 0.3])
    ax.tick_params(labelsize=18)
    plt.savefig(fout)

