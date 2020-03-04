voronoi_dict = None

gen_vmap_sigmap_params = dict(npix_per_reff=5,
                                  rscale=5.0,
                                  n_pseudo=1,
                                  verbose=False,
                                  voronoi=voronoi_dict,
                                  weight="luminosity",
                                  plot_map=True)

cal_lambda_params = dict(npix_per_reff=5,
                         rscale=5.0,
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
     'mstar_min': 5e8,
     'rmin': -1,
     'save': False,
     'unit_conversion': 'code',
     'verbose': False}

sfr_params = dict(hist_dt=0.1,
                  hist_tmin=0,
                  hist_tmax=None,
                  sfr_dts = [0.1, 0.5, 1.0])

gas_params = dict(dr=5, rmax=200, density_ratio=1e-3)
