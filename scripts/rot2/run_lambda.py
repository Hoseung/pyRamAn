import collections
from galaxymodule import galaxy
from galaxymodule import mk_gal_params as mgp
from galaxymodule import make_gal
from galaxymodule import vmax_sig
import pickle
import load
from load.hydro import Hydro
import tree
from analysis.cal_lambda import *
from galaxymodule.quick_mock import Simplemock
from general import defaults
from galaxymodule import rotation_parameter, rd_GM
import os
from utils.cosmology import Timeconvert

dfl = defaults.Default()

out_dir = "./gal_results_tree/"
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

# parameters
nout=782
Mcut = 1e10

gen_vmap_sigmap_params = dict(npix_per_reff=5,
                              rscale=3.0,
                              n_pseudo=60,
                              verbose=False,
                              voronoi=None, #voronoi_dict
                              weight="luminosity")


cal_lambda_params = dict(npix_per_reff=5,
                         rscale=3.0,
                         method='ellip',
                         verbose=False,
                         iterate_mge = False,
                         save_result = True,
                         galaxy_plot_dir='./')
mk_gal_params = dict()
mgp.HAGN["verbose"] = False


s = load.sim.Sim(nout=nout)
#gcat = tree.halomodule.Halo(nout=nout, is_gal=True)


#i_massive = np.where((gcat.data["id"] % 10 == 5) * (gcat.data["m"] > Mcut))[0]
#print("there are {} galaxies in the sample.".format(len(i_massive)))
#large_gals = gcat.data[i_massive]

MockSED = Simplemock(repo=dfl.dir_repo+"sed/")
timeconverter = Timeconvert(s.info)



#for large_gal in large_gals[100:200]:
if True:
    maintree, pidx, large_gals = pickle.load(open("all_direct_prgs_gal/ptrees/all_direct_prgs_7988258_gal.pickle", "rb"))
	
    gg = rd_GM.Gal(782, catalog=maintree[0], info=s.info)

    gg.debug=False
    #gg.mk_gal(**mgp.HAGN)
    make_gal.mk_gal(gg,**mgp.HAGN)
    gg.star['time'] = timeconverter.time2gyr(gg.star['time'],
                                 z_now = s.info.zred)
    gg.star.Flux_u = MockSED.get_flux(star=gg.star, filter_name="u")
    gg.star.Flux_r = MockSED.get_flux(star=gg.star, filter_name="r")
    # g-r color 
    # need correct unit and .... log scale subtract
    gg.meta.mean_age = np.mean(gg.star["time"])
    rotation_parameter.gen_vmap_sigmap(gg, **gen_vmap_sigmap_params)
    rotation_parameter.cal_lambda_r_eps(gg, **cal_lambda_params)

    # Calculate Vmax, Sig
    # get_vmax_sig uses vmap and sigmap from get_vmap_sigmap.
    # So if the maps were luminosity weighted, that also applies to the vmax and sig.
    vmax_sig.get_vmax_sig(gg, make_plot=False)

    # Save lambda results
    pickle.dump(gg.meta,
                open(out_dir + str(nout) + "_" + str(gg.meta.id) + \
                     str(gen_vmap_sigmap_params["n_pseudo"]) + \
                    ".pickle", "wb"))

