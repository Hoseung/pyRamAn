import pickle
import numpy as np
import matplotlib.pyplot as plt
#from galaxymodule.galaxy import Galaxy
import tree
import collections
#from galaxymodule import galaxy
import load
from galaxymodule import mk_gal_params as mgp
from analysis.cal_lambda import *
from make_gal import *

nout=782

nout_first = 310
Mcut = 1e10

# Load data

#s = load.sim.Sim(nout=nout)
#gcat = tree.halomodule.Halo(nout=nout, is_gal=True)
print("Start")

 
load_prg = True
dump_prg = False
 
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

nnza = np.genfromtxt("./nout_nstep_zred_aexp.txt",
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])

nnza["nstep"] = 62 - nnza["nstep"]
#nnza["nout"][:-1] = nnza["nout"][1:] # skip last snapshot

# Only for the later time

#
all_gcats=[]
all_infos=[]
nouts_info_cat=[]
nout_end = 733
for nout in nnza["nout"]:
    if nout > nout_first and nout <= nout_end:
        all_gcats.append(tree.halomodule.Halo(nout=nout, is_gal=True))
        all_infos.append(load.info.Info(nout=nout))
        nouts_info_cat.append(nout)

for fid,fidx in zip(final_ids, final_idxs):
    if load_prg:
        adp = pickle.load(open("./all_direct_prgs/all_direct_prgs_"+str(fidx)+".pickle", "rb"))
    else:
        maintree, idx_prgs_alltime = tt.extract_direct_full_tree(fidx)
        adp = tt.get_all_trees(idx_prgs_alltime)
        if dump_prg: 
            pickle.dump(adp, open("./all_direct_prgs/all_direct_prgs_"+str(fidx)+".pickle", "wb"))
        
    all_data=[]
    #main_gal=[]
    for sats_in_this_step in adp[1:]:
        # only for the mergers at last nout.
        # nout_now = nnza["nout"][np.abs(aexp - nnza["aexp"]).argmin()]
       
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
            for this_sat in sats_in_this_step:
                if this_sat is None:
                    continue
                print("Going through a satellite, # steps = ", len(this_sat))
                # progenitor of each satellite
                for istep, sat in enumerate(this_sat): 
                    nstep_sat = sat["nstep"]
                    inout_sat =np.where(nnza["nstep"] == nstep_sat)[0][0]
                    nout_sat = nnza["nout"][inout_sat]

                    print("SAT -progenitor \n", nstep_sat, nout_sat)
                    if nout_sat < nout_first:
                        continue
                    print(istep, inout)

                    #nout_sat = all_nouts[istep + inout] # from 782 to 43
                    #only at the last moment
                    rel_pos = (gg.header["xg"] - sat["xp"])*1e3 # in kpc
                    rel_vel = gg.header["vg"] - sat["vp"] 
                    jx,jy,jz = np.cross(rel_pos, rel_vel)
                    j_orbital=(jx,jy,jz)/np.sqrt(jx**2 + jy**2 + jz**2)
         
                    # spin alignment
                    gcat_now = all_gcats[i_info_cat + istep]
                    info_now = all_infos[i_info_cat + istep]
                    try:
                        print(sat["xp"])
                        gsat = load.rd_GM.Gal(nout=nout_sat,
                               catalog=np.copy(gcat_now.data[sat["id"]-1]),
                               info=info_now)    
                        gsat.debug = False
                        mk_gal(gsat, **mgp.HAGN)
         
                        gsat.cal_norm_vec()
                        print(gsat.meta.xc)
                        # orbit
                        gsat.meta.j_orbit = (jx,jy,jz)
                        gsat.meta.relang = 180 / np.pi * np.arccos(np.dot(gg.meta.nvec, j_orbital))
                        # spin
                        gsat.meta.spinang = 180./np.pi*np.arccos(np.dot(gg.meta.nvec, gsat.meta.nvec))
                    except:
                        print("failed : ",sat["id"])
                        #this_sat_spinang.append(-1)
                        gsat.meta.spinang = -1.
                    # tree information
                    gsat.meta.root_id = fid
                    gsat.meta.root_idx = fidx
                    gsat.meta.root_sat_id = this_sat["id"]
                    gsat.meta.root_sat_idx = this_sat["idx"]
                    gsat.meta.nstep_merge = nnza["nstep"][inout]
                    gsat.meta.nout_merge = nout_now
                    gsat.meta.nstep_sat_merge = nnza["nstep"][istep + inout]
                    gsat.meta.nout_sat_merge = nout_sat
                    all_data.append(gsat.meta) 
    #all_data.extend(main_gal) 
    pickle.dump(all_data, open("all_meta_data"+str(fidx)+".pickle", "wb"))
