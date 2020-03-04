import numpy as np

def add_output_containers(gg):
    gg.meta.sfr_results={"hist_dt":None, "hist_tmin":None, "hist_tmax":None, "hist":None, "sfr_dts":None, "sfrs":None, "area":None}
    #gg.meta.lambda_results={"lambda_results", }
    #gg.meta.mge_results={"mge_results":None}
    gg.meta.gas_results={"gas_results":None, "mgas_tot":None, "mgas_cold":None, "Ln_gas":None}
    gg.meta.vsig_results={"Vmax":None, "sigma":None, "V_sig":None}

def ind_cell_kd(kdtree, gal, pboxsize, rscale=25.0, rmax=500, rmin=30):
    """
    Extract cells within rscale * Rreff and add to the galaxy.
    """
    xc,yc,zc = gal.meta.xc, gal.meta.yc, gal.meta.zc
    rgal = min([rmax, max([rmin, gal.meta.reff * rscale])]) / (pboxsize*1e3) # kpc -> code unit
    #index = kdtree.query_ball_point((xc,yc,zc), rgal)
    xyzcen = (xc/pboxsize + 0.5,
              yc/pboxsize + 0.5,
              zc/pboxsize + 0.5)
    return kdtree.query_ball_point(xyzcen, rgal)


def get_cell(allcell, kdtree, gg, info, **kwargs):
    # Simple spherical cut.
    gg.cell=allcell[ind_cell_kd(kdtree, gg, info.pboxsize, **kwargs)]

    if len(gg.cell) > 1:
        #print(s.hydro.cell["x"].ptp())
        # convert to kpc
        #print("gg.center_code", gg.center_code)
        gg.cell["pos"] = (gg.cell["pos"]-gg.center_code)*info.boxtokpc
        gg.cell["vel"] = gg.cell["vel"]*info.kms -gg.meta.vel
        gg.cell["dx"] *= info.boxtokpc
