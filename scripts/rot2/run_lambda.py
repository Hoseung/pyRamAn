from galaxymodule import galaxy
from galaxymodule import mk_gal_params as mgp
from galaxymodule import make_gal
from galaxymodule import vmax_sig
import pickle
import load
from load.hydro import Hydro
import tree
#from analysis.cal_lambda import *
from galaxymodule.quick_mock import Simplemock
from general import defaults
from galaxymodule import rotation_parameter, rd_GM
import os
from utils.cosmology import Timeconvert
from utils import hagn
from glob import glob
import numpy as np



dfl = defaults.Default()
nnza = hagn.Nnza()


out_dir = "./gal_results_tree/"
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

# parameters
nout=nnza.nnza["nout"].max()
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

#s = load.sim.Sim(nout=nout)
#gcat = tree.halomodule.Halo(nout=nout, is_gal=True)


#i_massive = np.where((gcat.data["id"] % 10 == 5) * (gcat.data["m"] > Mcut))[0]
#print("there are {} galaxies in the sample.".format(len(i_massive)))
#large_gals = gcat.data[i_massive]

# All infos
all_infos=[]
for nout in nnza.nnza["nout"]:
    all_infos.append(load.info.Info(nout=nout))



MockSED = Simplemock(repo=dfl.dir_repo+"sed/")


timeconverter = Timeconvert(all_infos[0])


all_files = glob("all_direct_prgs_gal/ptrees/all_direct*.pickle")
#for large_gal in large_gals[100:200]:
for this_file in all_files[600:800]:
    maintree, pidx, sats = pickle.load(open(this_file, "rb"))
    result_gal_evol=[]
    for istep, this_main in enumerate(maintree):
        gid = this_main["id"]
        if gid <= 0:
            continue
        mgp.HAGN["mstar_min"]=1e9
        info = all_infos[istep]
        nout = nnza.step2out(this_main["nstep"])
        # GalaxyMaker files are not available....!!
        if nout < 303:
            continue
        gg = rd_GM.Gal(nout,
                      catalog=this_main.copy(), info=info)
 
        gg.debug=False
        #gg.mk_gal(**mgp.HAGN)
        make_gal.mk_gal(gg,**mgp.HAGN)
        gg.star['time'] = timeconverter.time2gyr(gg.star['time'],
                                     z_now = info.zred)
        #print(gg.star["time"].min(), gg.star["time"].max())
        gg.star.Flux_u = MockSED.get_flux(star=gg.star, filter_name="u")
        gg.star.Flux_r = MockSED.get_flux(star=gg.star, filter_name="r")
        # g-r color 
        # need correct unit and .... log scale subtract
        gg.meta.mean_age = np.mean(gg.star["time"])
        gg.cal_norm_vec()

        #Below this line # part is high.
        rotation_parameter.gen_vmap_sigmap(gg, **gen_vmap_sigmap_params)
        rotation_parameter.cal_lambda_r_eps(gg, **cal_lambda_params)
 
        # Calculate Vmax, Sig
        # get_vmax_sig uses vmap and sigmap from get_vmap_sigmap.
        # So if the maps were luminosity weighted, that also applies to the vmax and sig.
        vmax_sig.get_vmax_sig(gg, make_plot=False)


        gg.meta.nout = nout
        gg.meta.nstep = this_main["nstep"]
        gg.meta.idx = this_main["idx"]
 
        result_gal_evol.append(gg.meta)

    # Save lambda results
    pickle.dump(result_gal_evol,
                open(out_dir + str(maintree[0]["idx"]) + \
                str(gen_vmap_sigmap_params["n_pseudo"]) + \
                        ".pickle", "wb"))

