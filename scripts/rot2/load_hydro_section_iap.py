from cal_gas import *
import load
import pickle
import numpy as np
import tree
from makegal import *
from galaxymodule import mk_gal_params as mgp

nout = 782
s = load.sim.Sim(nout=nout, setup=True)
gcat = tree.halomodule.Halo(nout=nout, is_gal=True)

nnza = np.genfromtxt("./nout_nstep_zred_aexp.txt",
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])


all_direct_prgs = pickle.load(open("all_direct_prgs_7864636.pickle", "rb"))
ids, pos, rvir, nstep = get_flat_data(all_direct_prgs, ["id", "xp","rvir","nstep"])
nstep_now = 61
pos_now = pos[nstep == nstep_now]
rvir_now = rvir[nstep == nstep_now]
ids_now = ids[nstep == nstep_now]


ranges = range_from_pos_r(pos_now, rvir_now, rscale=1.0, buffer=0.1)
ranges = range_2code_unit(ranges, s.info.pboxsize)
print(ranges)
s.set_ranges(ranges)

s.add_hydro()


###
# mk_gal using cell data 
# Cell to be in the GalaxyMaker coordinate.
# 
#
cell = s.hydro.cell
print(len(cell))

nout_now = nnza["nout"][np.where(61 - nnza["nstep"] == nstep_now)[0]][0]
for gid in ids_now:
    gg = load.rd_GM.Gal(nout=nout_now,
                    catalog=np.copy(gcat.data[gid-1]),
                    info=s.info)
    #print("gg.star org", gg.star["x"])
    gg.cell = cell # Will it be copied internally? 
    #gg.__has_cell__=True
    gg.debug = False
    mk_gal(gg, **mgp.HAGN)
    gg.cal_norm_vec()
    gg.cal_mgas()
###


