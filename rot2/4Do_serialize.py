import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import galaxymodule  # needed for result_sub_sample_**.pickle
import utils.match as mtc
import pickle
#from glob import glob

from utils import hagn
import os
from rot2.analysis import *
from rot2 import serialize_results
import tree.halomodule as hmo
from rot2 import cell_chunk_module as ccm

from rot2 import serialize_results
from utils import cosmology
from load import info



out_base="./RUN1/"

prg_dir = out_base+"test_direct_prgs_gal/"
result_dir = out_base+"lambda_results/"

# Load prg tree data

all_ids= np.genfromtxt(prg_dir + "final_idxs_allmassive_gal.txt", dtype=int)
all_final_idxs = all_ids[:,0]
all_final_ids = all_ids[:,1]

nnza_all = hagn.Nnza(fname=out_base+"nout_nstep_zred_aexp.txt")
nnza_cell = hagn.Nnza(fname=out_base+"nout_nstep_zred_aexp_63.txt")

istep_max = 40
nouts = nnza_cell.nnza["nout"][:istep_max]

allresults = get_all_results(nouts, prg_dir=prg_dir, out_dir =result_dir, fix_idx_nout=9654)



all_sample_idxs=pickle.load(open(prg_dir + "all_sample_idxs.pickle", "rb"))


serial_out_dir = prg_dir+"result_serial/"
if not os.path.isdir(serial_out_dir):
    os.mkdir(serial_out_dir)

# Build serial results and dump.
# Chunks of results in each nout (lambda_results/{nout}/result_sub_sample_{nout})

Allallidxs=[]
for result_thisnout in allresults:
    Allallidxs.append(np.array([agal.idx for agal in result_thisnout]))
Allallids=[]
for result_thisnout in allresults:
    Allallids.append(np.array([agal.id for agal in result_thisnout]))

info = info.Info(nout=782)

tc = cosmology.Timeconvert(info, zred_now=0)

import numpy.lib.recfunctions as recf

all_fid_ok=[]
all_fidx_ok=[]

nstep_too_short_main = 100
for i, this_idx in enumerate(all_final_idxs):
    fname = prg_dir + "{}_adp.pickle".format(this_idx)
    if not os.path.isfile(fname):
        # dump_prgs    broken in the middle.
        continue
    #all_fid_ok.append(this_id)
    #print(i, "IDX=",this_idx)
    adp = pickle.load(open(fname, "rb"))
    if min(adp[0][0]["nstep"][adp[0][0]["nstep"] > 0]) > nstep_too_short_main:
        print("skip")
        continue
    #print(len(adp))
    #print(hasattr(this_gal.maintree,"mask"))
    # Append age to maintree and mainresult.
    lbt = tc.zred2gyr(nnza_all.a2b(adp[0][0]["nstep"],"nstep","zred"),z_now=0)
    adp[0][0] = recf.append_fields(adp[0][0], "time", lbt)
    max_step=len(allresults)
    this_gal = serialize_results.Serial_result(adp)
    cnt_merger=0
    for i, this_sats in enumerate(adp):
        nout_results=[]
        #print("# sats", len(this_sats))
        for sat in this_sats:
            sat_results=[]

            #if len(sat) > 0 :
                #print("merger1")
            for ss in sat:
                nout=nnza_all.step2out([ss["nstep"]]) # NOT nnza_cell.
                if nout in nouts:
                    istep_cell = np.where(nnza_cell.nnza["nout"] == nout)[0][0]
                    allresults_now=allresults[istep_cell]
                    allresults_now_idx=Allallidxs[istep_cell]
                    i_result = np.where(allresults_now_idx == ss["idx"])[0]
                    #print("len results", len(allresults_now), "i_result", i_result)
                    if len(i_result) > 0:
                        sat_results.append(allresults_now[i_result[0]])
                        sat_results[-1].nout=int(nout)
                        sat_results[-1].nstep=nnza_cell.nnza["nstep"][istep_cell]
            #print(len(sat_results))
            nout_results.append(sat_results)
            # Merger properties
            if i == 0:
                this_gal.main_arr = serialize_results.galresult2rec(sat_results)
            elif len(sat_results) > 0 and sat_results[0].mstar > 0.0:
                #print("merger2")
                this_gal.add_merger(sat_results, sat)
                cnt_merger+=1
                #this_gal.mergers.append(serialize_results.get_merger_props(this_gal.main_arr,
                #                        serialize_results.galresult2rec(sat_results)))
        #print(i)
        this_gal.data.append(nout_results)


    pickle.dump(this_gal, open(serial_out_dir+"serial_result{}.pickle".format(this_idx), "wb"))
    all_fidx_ok.append(this_idx)
    #break


np.savetxt(serial_out_dir+"all_fidx_ok.txt", all_fidx_ok, fmt="%d")
serial_results = [pickle.load(open(serial_out_dir+"serial_result{}.pickle".format(this_idx), "rb")) for this_idx in all_fidx_ok]
