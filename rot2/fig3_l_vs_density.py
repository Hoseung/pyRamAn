from analysis.all_plot_modules import density_map
import numpy as np
import matplotlib.pyplot as plt

def lambda_vs_density(allresults, out_dir="./"):

    den_cmap = "YlGnBu_r"
    allids = [gal.id for gal in allresults]

    allids=[]
    lambda_r=[]
    eps=[]
    mgal=[]
    mgas_c=[]
    mgas_t=[]
    very_young_stars=[]
    young_stars=[]
    for gal in allresults:
        if hasattr(gal, "mge_result_list"):
            eps.append(gal.mge_result_list[0]["eps"])
            try:
                lambda_r.append(gal.lambda_r[0])
            except:
                lambda_r.append(gal.lambda_r)
            allids.append(gal.id)

            mgal.append(gal.mstar)
            mgas_c.append(gal.gas_results["mgas_cold"])
            mgas_t.append(gal.gas_results["mgas_tot"])
            young_stars.append(sum(gal.sfr_results["hist"][:10])*gal.sfr_results["area"] * 1e9 * gal.sfr_results["hist_dt"])
            very_young_stars.append(gal.sfr_results["hist"][0]*gal.sfr_results["area"] * 1e9 * gal.sfr_results["hist_dt"])

    allids = np.array(allids)
    eps = np.array(eps)
    lambda_r = np.array(lambda_r)
    mgal = np.array(mgal)
    mgas_c = np.array(mgas_c)
    mgas_t = np.array(mgas_t)
    young_stars = np.array(young_stars)
    vert_young_stars = np.array(very_young_stars)

    fig, axs = plt.subplots(1,3, sharey=True)

    #axs[0].scatter(np.log10(mgas_t/mgal+1), lambda_r, alpha=0.3)
    #axs[1].scatter(np.log10(young_stars/mgal+1), lambda_r, alpha=0.3)
    #axs[2].scatter(np.log10(very_young_stars/mgal+1), lambda_r, alpha=0.3)
    #plt.savefig()

    xx,yy,z = density_map(np.log10(mgas_c/mgal+1), lambda_r)
    axs[0].scatter(xx, yy, c=z, s=30, edgecolor='',
                   cmap=den_cmap, rasterized=False,
                   alpha=1.0)
    xx,yy,z = density_map(np.log10(mgas_t/mgal+1), lambda_r)
    axs[1].scatter(xx, yy, c=z, s=30, edgecolor='',
                   cmap=den_cmap, rasterized=False,
                   alpha=1.0)

    xx,yy,z = density_map(np.log10(mhost), lambda_r)
    axs[2].scatter(xx, yy, c=z, s=30, edgecolor='',
                   cmap=den_cmap, rasterized=False,
                   alpha=1.0)

    axs[0].set_ylabel(r"$\lambda_{Re}$")
    axs[0].set_xlabel("D10 (>1e10) in pMpc")
    axs[1].set_xlabel("D50 (>1e10) in pMpc")
    axs[2].set_xlabel(r"log($M_{host}$)")

    plt.tight_layout()
    plt.savefig(out_dir + "lambda_vs_density.png", dpi=200)
