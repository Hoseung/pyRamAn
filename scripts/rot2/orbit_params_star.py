import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import Arrow
#from galaxymodule.galaxy import Galaxy
import tree
import collections
#from galaxymodule import galaxy
import load
from galaxymodule import mk_gal_params as mgp
from analysis.cal_lambda import *
from makegal import *

def all_cats_infos(nouts, is_gal=True):
    all_gcats=[]
    all_infos=[]
    nouts_info_cat=[]
    for nout in nouts:
        all_gcats.append(tree.halomodule.Halo(nout=nout, is_gal=True))
        all_infos.append(load.info.Info(nout=nout))
        nouts_info_cat.append(nout)
    return all_gcats, all_infos, nouts_info_cat

class Nnza():
    def __init__(self, fname="./nout_nstep_zred_aexp.txt"):
        self.nnza = np.genfromtxt(fname,
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])
        #self.nnza["nstep"] = 62 - self.nnza["nstep"]
    def step2out(self, nstep):
        return int(self.nnza["nout"][np.where(self.nnza["nstep"] == nstep)[0]])

    def out2step(self, nout):
        return int(self.nnza["nstep"][np.where(self.nnza["nout"] == nout)[0]])


nout=782
nout_first = 400

# Load data
# Draw stellar density map of host and a satellite.
plot_stars = False
# Load existing progenitor trees. 
# Much faster by skipping reading tree data.
load_prg = True
# save progenitor trees 
dump_prg = False
 
nnza = Nnza(fname="./nout_nstep_zred_aexp.txt")
nnza = nnza.nnza


if not load_prg:
    tt = tree.tmtree.Tree(is_gal=True)
    print("Loading tree done")

    tnow = tt.tree[tt.tree["nstep"]==max(tt.tree["nstep"])]
    large_last = tnow[(tnow["m"] > 1) * (tnow["m"] < 2)]

    final_idxs = large_last["idx"][large_last["id"] % 10 == 5]
    final_ids = large_last["id"][large_last["id"] % 10 == 5]

    np.savetxt("final_idxs.txt", final_idxs, fmt='%d')
    np.savetxt("final_ids.txt", final_ids, fmt='%d')
else: 
    final_idxs = np.genfromtxt("final_idxs.txt", dtype="<i4")
    final_ids = np.genfromtxt("final_ids.txt", dtype="<i4")

# Needed to convert between nout and nstep.
# gcat is needed for reading an individual galaxy file.
# Because the code traverses along each satellite tree,
# catalogs at all nouts are accessed repeatedly.
nout_end = 733

    
    

#for nout in nnza["nout"]:
#    if nout > nout_first and nout <= nout_end:
all_gcats, all_infos, nouts_info_cat = all_cats_infos(\
            nnza["nout"][(nnza["nout"] > nout_first) * (nnza["nout"] <= nout_end)])


# Main loop
for fid,fidx in zip(final_ids, final_idxs):
    if load_prg:
        adp = pickle.load(open("./all_direct_prgs/all_direct_prgs_"+str(fidx)+".pickle", "rb"))
    else:
        maintree, idx_prgs_alltime = tt.extract_direct_full_tree(fidx)
        adp = tt.get_all_trees(idx_prgs_alltime)
        if dump_prg: 
            pickle.dump(adp, open("./all_direct_prgs/all_direct_prgs_"+str(fidx)+".pickle", "wb"))
        
    all_data=[]

    for sats_in_this_step in adp:
        # adp[0] = all mergers at the last snapshot.
        # adp[1] = all mergers at the last -1 snapshot.
        if len(sats_in_this_step) == 0:
            continue
        if sats_in_this_step[0] == None:
            continue
        nstep_this_merge = sats_in_this_step[0][0]["nstep"]
        id_host_this_step = sats_in_this_step[0][0]["id"]
        inout = np.where(nnza["nstep"] == nstep_this_merge)[0][0]
        nout_now=nnza["nout"][inout]
        if nout_now < nout_first:
            continue
        i_info_cat = int(np.where(nouts_info_cat == nout_now)[0])
        gcat = all_gcats[i_info_cat]
        gg = load.rd_GM.Gal(nout=nout_now,
                            catalog=np.copy(gcat.data[id_host_this_step-1]),
                            info=all_infos[i_info_cat])
        gg.debug=False
        mgp.HAGN["mstar_min"] = 1e9
        mk_gal(gg,**mgp.HAGN)
        gg.cal_norm_vec()
        # Make sure that the catalog is not modified.
        
        gg.meta.root_id = fid
        gg.meta.root_idx = fidx
        gg.meta.root_sat_id = -1
        gg.meta.root_sat_idx = -1
        gg.meta.nstep_merge = -1
        gg.meta.nout_merge = -1
        gg.meta.nstep_sat_merge = -1
        gg.meta.nout_sat_merge = -1
        all_data.append(gg.meta)
        if len(sats_in_this_step) > 1:
             # Has merger (0 = self)
            if plot_stars:
                fig, ax = plt.subplots()
            for this_sat in sats_in_this_step:
                if this_sat is None:
                    continue
                print("Going through a satellite", len(this_sat))
                # progenitor of each satellite
                for istep, sat in enumerate(this_sat):
                    print("In a satellite")
                    # only at the last moment
                    nstep_sat = sat["nstep"]
                    inout_sat =np.where(nnza["nstep"] == nstep_sat)[0][0]
                    nout_sat = nnza["nout"][inout_sat]

                    if nout_sat < nout_first:
                        continue

                    #nout_sat = all_nouts[istep + inout] # from 782 to 43
                    #only at the last moment
                    rel_pos = (gg.header["xg"] - sat["xp"])*1e3 # in kpc
                    rel_vel = gg.header["vg"] - sat["vp"] 
                    jx,jy,jz = np.cross(rel_pos, rel_vel)
                    #print("JXJYJZ", jx,jy,jz)
                    j_orbital=(jx,jy,jz)/np.sqrt(jx**2 + jy**2 + jz**2)
         
                    # spin alignment
                    gcat_now = all_gcats[i_info_cat + istep]
                    info_now = all_infos[i_info_cat + istep]
                    if True:
                        #print(sat["xp"])
                        gsat = load.rd_GM.Gal(nout=nout_sat,
                               catalog=np.copy(gcat_now.data[sat["id"]-1]),
                               info=info_now)    
                        gsat.meta.j_orbit = (jx,jy,jz)
                        gsat.meta.relang = 180 / np.pi * np.arccos(np.dot(gg.meta.nvec, j_orbital))
                        gsat.debug = False
                        mk_gal(gsat, **mgp.HAGN)
         
                        gsat.cal_norm_vec()
                        #print(gsat.meta.xc)
                        # orbit
                        # spin
                        print("gsat done")
                        gsat.meta.spinang = 180./np.pi*np.arccos(np.dot(gg.meta.nvec, gsat.meta.nvec))
                    if False:
                        print("failed : ",sat["id"])
                        #this_sat_spinang.append(-1)
                        gsat.meta.spinang = -1.
                    # tree information
                    gsat.meta.root_id = fid
                    gsat.meta.root_idx = fidx
                    gsat.meta.root_sat_id = this_sat["id"][0]
                    gsat.meta.root_sat_idx = this_sat["idx"][0]
                    gsat.meta.nstep_merge = nnza["nstep"][inout]
                    gsat.meta.nout_merge = nout_now
                    gsat.meta.nstep_sat_merge = nnza["nstep"][istep + inout]
                    gsat.meta.nout_sat_merge = nout_sat
                    all_data.append(gsat.meta) 
                    
                    if plot_stars:
                        xc_pri = gg.header["xg"][0]*1e3
                        yc_pri = gg.header["xg"][1]*1e3
                        xc_sec = gsat.header["xg"][0]*1e3
                        yc_sec = gsat.header["xg"][1]*1e3
                        #print("Before map", gg.star["x"])
                        print("center", xc_pri)
                        ax.hist2d(np.concatenate((gg.star["x"] + xc_pri,
                                                  gsat.star["x"] + xc_sec)),
                                  np.concatenate((gg.star["y"] + yc_pri,
                                                  gsat.star["y"] + yc_sec)),
                                  weights=np.concatenate((gg.star["m"], gsat.star["m"])),
                                  bins=200, range=[[xc_pri-5e2, xc_pri+5e2],
                                                   [yc_pri-5e2, yc_pri+5e2]])
                        
                        # spin vector of the primary
                        ax.add_patch(Arrow(xc_pri,yc_pri,
                                           gg.meta.nvec[0]*5e1,
                                           gg.meta.nvec[1]*5e1, width=4.0 ))
                        # spin vector of the secondary
                        ax.add_patch(Arrow(xc_sec,yc_sec,
                                     gsat.meta.nvec[0]*5e1,
                                     gsat.meta.nvec[1]*5e1, width=4.0)) 
                        ax.set_aspect('equal', 'datalim')
                        plt.savefig("gals_2d_"+str(fidx) +"_"+str(gsat.meta.root_sat_idx)+"_"+
                                    str(gsat.meta.nstep_sat_merge) + ".png")
                        plt.cla()


    #all_data.extend(main_gal) 
    pickle.dump(all_data, open("all_meta_data"+str(fidx)+".pickle", "wb"))
