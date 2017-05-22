from cal_gas import *
import load
import pickle
import tree
from makegal import *
from galaxymodule import mk_gal_params as mgp


nout = 782
s = load.sim.Sim(nout=nout, setup=True)
gcat = tree.halomodule.Halo(nout=nout, is_gal=True)


if 1 == 2:
    all_direct_prgs = pickle.load(open("all_direct_prgs_7864636.pickle", "rb"))
    ids, pos, rvir, nstep = get_flat_data(all_direct_prgs, ["id", "xp","rvir","nstep"])
    nstep_now = 61
    pos_now = pos[nstep == nstep_now]
    rvir_now = rvir[nstep == nstep_now]
else:
    nstep_now = 62
    pos_now = gcat.data["x"]
    rvir_now = gcat.data["x"]

#ranges = range_from_pos_r(pos_now, rvir_now, rscale=1.5, buffer=0.5)

gid =125
xc = gcat.data[gid - 1]["x"]
yc = gcat.data[gid - 1]["y"]
zc = gcat.data[gid - 1]["z"]
r =  gcat.data[gid - 1]["rvir"]
rs = 2
#print(xc,yc,zc,r)
ranges = [[xc-r*rs,xc+r*rs],
          [yc-r*rs,yc+r*rs],
          [zc-r*rs,zc+r*rs]]

#print(ranges)
s.set_ranges(ranges)
#print(s.ranges)
s.add_hydro()

print("len(cell)", len(s.hydro.cell))
###
# mk_gal using cell data 
# Cell to be in the GalaxyMaker coordinate.
#  
#
cell = s.hydro.cell
print(cell["x"])
if 2 == 2:
    cell["x"] -= xc
    cell["y"] -= yc
    cell["z"] -= zc
    cell["x"] *= s.info.boxtokpc
    cell["y"] *= s.info.boxtokpc
    cell["z"] *= s.info.boxtokpc
    cell["dx"] *= s.info.boxtokpc

# O, Fe, C, N, Mg, Si are tracked!

#for gid in ids:
gg = load.rd_GM.Gal(nout=nout,
                   catalog=np.copy(gcat.data[gid-1]),
                   info=s.info)
gg.cell = cell # Will it be copied internally? 
print(gg.star["x"])
print(gg.cell["x"])
gg.debug = False
mgp.HAGN["verbose"] = True
mk_gal(gg, **mgp.HAGN)
gg.cal_norm_vec()
    
###


