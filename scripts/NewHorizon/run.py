import matplotlib
matplotlib.use("qt5agg")

import analysis.NH_module as nhm
import tree.halomodule as hmo
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from galaxymodule import make_gal
from galaxymodule import mk_gal_params as mgp
import utils

from utils.sampling import Region
from utils.cosmology import Timeconvert
from galaxymodule import rotation_parameter, rd_GM
from scipy.spatial import cKDTree
from load.sim import Sim
from galaxymodule import vmax_sig
import pickle
import galaxymodule.quick_mock as qmc
from galaxymodule import gal_properties
import pickle
import os
from galaxymodule import rotation_parameter as rotp


# Load a galaxy

# Parameters----------------------------------------------
nout=606
out_base='./'
wdir = "./"
rscale=1.5
do_cell=True # Needed to measure the total mass
do_color=False
out_subdir = "lambda_results/"
save_cell=True
do_plot=True
voronoi_dict = None

gen_vmap_sigmap_params = dict(npix_per_reff=15,
                                  rscale=3.0,
                                  n_pseudo=1,
                                  verbose=False,
                                  voronoi=voronoi_dict,
                                  weight="mass",
                                  plot_map=True)

cal_lambda_params = dict(npix_per_reff=15,
                         rscale=3.0,
                         method='ellip',
                         verbose=False,
                         voronoi=voronoi_dict,
                         iterate_mge = False,
                         save_result = True,
                         galaxy_plot_dir='./',
                         recenter_v=True)

mgp_NH = {'Rgal_to_reff': 5.0,
     'den_lim': 1e5,
     'den_lim2': 3.333e5,
     'follow_bp': None,
     'method_com': 1,
     'mstar_min': 1e9,
     'rmin': -1,
     'save': False,
     'unit_conversion': 'code',
     'verbose': False}


# Basic setup---------------------------------------------
out_dir = out_base+out_subdir+ str(nout) +'/'


# Load a galaxy

# Parameters----------------------------------------------
nout=606
out_base='./'
wdir = "./"
rscale=1.5
do_cell=True # Needed to measure the total mass
do_color=False
out_subdir = "lambda_results/"
save_cell=True
do_plot=True
voronoi_dict = None

gen_vmap_sigmap_params = dict(npix_per_reff=15,
                                  rscale=3.0,
                                  n_pseudo=1,
                                  verbose=False,
                                  voronoi=voronoi_dict,
                                  weight="mass",
                                  plot_map=True)

cal_lambda_params = dict(npix_per_reff=15,
                         rscale=3.0,
                         method='ellip',
                         verbose=False,
                         voronoi=voronoi_dict,
                         iterate_mge = False,
                         save_result = True,
                         galaxy_plot_dir='./',
                         recenter_v=True)

mgp_NH = {'Rgal_to_reff': 5.0,
     'den_lim': 1e5,
     'den_lim2': 3.333e5,
     'follow_bp': None,
     'method_com': 1,
     'mstar_min': 1e9,
     'rmin': -1,
     'save': False,
     'unit_conversion': 'code',
     'verbose': False}


# Basic setup---------------------------------------------
out_dir = out_base+out_subdir+ str(nout) +'/'


import analysis.gal_props as gp
import galaxymodule
import importlib
importlib.reload(gp)
importlib.reload(nhm)
importlib.reload(galaxymodule)
importlib.reload(galaxymodule.rd_GM)
importlib.reload(galaxymodule.galaxy)


do_dm = True
do_cell=True

band = "flux_u"

msun_in_g = 1.989e33
kpc_in_cm = 3.086e+21

new_dtype = {"ellip" : (('<f8', 1), 0),
             "flux_u": (('<f8', 1), 8),
             "flux_g": (('<f8', 1), 16),
             "flux_r": (('<f8', 1), 24),
             "flux_i": (('<f8', 1), 32)}

mgp.NH["method_cov"] = "close_member"
gas_params = dict(dr=5, rmax=200, density_ratio=1e-3)


#gid = 15
for gcat_this in good_gals[2:3]:
    gid = gcat_this["id"]
    print("ID gal =", gid)
    hid = gcat_this["hosthalo"]
    thishalo = hcat.data[hid-1].copy()

    gcat_sub = gcat.data[gcat.data["host"]==gid]
    hcat_sub = hcat.data[hcat.data["host"]==hid]
    gg = rd_GM.Gal(nout=nout, wdir=wdir,
                   gcat=gcat_this.copy(),
                   hcat=thishalo,
                   #gcat_subs = gcat_sub,
                   hcat_subs = hcat_sub,
                   info=s.info, type_cell=do_cell,
                   load=False)

    cell_to_msun = gg.info.unit_d/msun_in_g * kpc_in_cm**3

    gg.debug=False
    gg.meta.nout = nout
    gg.meta.idx = gcat_this["idx"]

    rd_star_params = dict(additional_fields=new_dtype)
    gg.load(type_cell="none", type_dm="none", rd_star_params=rd_star_params)

    if not make_gal.mk_gal(gg,**mgp_NH):
        print("Too small")

    #print(i,"-th, idx=", gg.meta.idx)
    if gg.meta.nstar < 500:
        gg.meta.mgas_cold = -1
        result_sub_sample.append(gg.meta)


    # Add output dicts
    nhm.add_output_containers(gg)

    nhm.load_components(gg, s,
                        load_dm=do_dm,
                        load_cell=do_cell,
                        load_raw=True,
                        save_cell=True,
                        idlist=hcat.idlists[thishalo["id"]-1])
    # Now star and cell memberships are determined.

    if do_cell:
        gal_properties.get_cold_cell(gg, s.info, **gas_params)
        gal_properties.get_gas_properties(gg, s.info)
    else:
        gg.cell = None


    gg.cal_norm_vec(pop_nvec=["star"])
    if do_cell:
        gg.reorient(dest=[0,1,0], pops=["star", "cell"])
    else:
        gg.reorient(dest=[0,1,0], pops=["star"])

    if do_cell:
        gmass = gg.cell["var0"] * gg.cell["dx"]**3 * kpc_to_cm**3 * s.info.unit_d / (2*1e33) # in Msun
        angh = (gg.cell["x"]*gg.cell["vy"] -gg.cell["y"] *gg.cell["vx"]) * gmass
        gg.meta.ang_h = np.sum(angh)

    angs = (gg.star["x"]*gg.star["vy"]-gg.star["y"] *gg.star["vx"])*gg.star["m"]
    gg.meta.ang_s = np.sum(angs)

    nhm.plot_2d_simple_maps(gg)
    if do_cell: nhm.plot_2d_simple_gas_maps(gg)
    nhm.plot_radial(gg)
    plt.show()

    # Circularity
    gp.get_E(gg, nvec=[0,1,0], phi_direct=True)
    # This is totally wrong!
    gg.cal_b2t(ptype="star", bound_only=False)
    #print(gg.meta.b2t)

    # Stellar luminosity
    if gg.star['time'].max() < 0:
        gg.star['time'] = tc.time2gyr(gg.star['time'], z_now = gg.info.zred)
    gg.star["flux_u"] = MockSED.get_flux(gg.star, cell=None, quick=True, speed_check=True, filter_name="u")
    gg.star["flux_g"] = MockSED.get_flux(gg.star, cell=None, quick=True, speed_check=True, filter_name="g")
    gg.star["flux_r"] = MockSED.get_flux(gg.star, cell=None, quick=True, speed_check=True, filter_name="r")
    gg.star["flux_i"] = MockSED.get_flux(gg.star, cell=None, quick=True, speed_check=True, filter_name="i")

    bulge = gg.star[gg.star["ellip"] < 0.5]
    disk = gg.star[gg.star["ellip"] > 0.8]
    nhm.plot_bulge_disk(gg, bulge, disk, band="flux_u")

    print("B/T (0.5) = {:.2f}".format(np.sum(bulge["m"])/np.sum(gg.star["m"])))
