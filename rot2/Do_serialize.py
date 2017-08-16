import numpy as np
import galaxymodule  # needed for result_sub_sample_**.pickle
import utils.match as mtc
import pickle
#from utils import hagn
import os
from rot2.analysis import *
from rot2 import serialize_results
#import tree.halomodule as hmo
#from rot2 import cell_chunk_module as ccm
import numpy.lib.recfunctions as recf
from utils import cosmology
from load import info

def fill_main(mainarr, nnza_cell):
    # Set up a new array.
    new_nouts = nnza_cell.nnza["nout"][:mainarr["nstep"].ptp()+1]
    newarr = np.zeros(len(new_nouts), dtype=mainarr.dtype)
    # It's easy to fill nouts and nsteps.
    newarr["nout"]=new_nouts
    newarr["nstep"]=nnza_cell.a2b(newarr["nout"], "nout", "nstep")

    interp_fields = list(this_gal.main_arr.dtype.names)
    for field in ["nout", "nstep"]:
        interp_fields.remove(field)

    lbt_org = tc.zred2gyr(nnza_cell.a2b(mainarr["nstep"],"nstep","zred"), z_now=0)
    lbt_new = tc.zred2gyr(nnza_cell.a2b(newarr["nstep"],"nstep","zred"), z_now=0)

    for field in ["id", "idx"]:
        newarr[field][mtc.match_list_ind(newarr["nout"], mainarr["nout"])] = mainarr[field]
        interp_fields.remove(field)

    for field in interp_fields:
        if mainarr[field].ndim == 2:
            for i in range(3):
                r_p = mainarr[field][:,i]
                newarr[field][:,i] = np.interp(lbt_new, lbt_org, mainarr[field][:,i])
        else:
            r_p = mainarr[field]
            newarr[field] = np.interp(lbt_new, lbt_org, r_p)
    return newarr

# interpolate main galaxy results on finetree.


def interpol_fine(this_gal, nnza_cell, nnza_all, do_smooth=True):
    finetree=this_gal.maintree
    mainarr = this_gal.main_arr
    finearr = np.zeros(len(finetree),dtype=mainarr.dtype)
    fields_interp = list(mainarr.dtype.names)

    finearr["nstep"]=finetree["nstep"]
    finearr["id"] = finetree["id"]
    finearr["idx"] = finetree["idx"]
    finearr["pos"] = finetree["pos"] # Pos and vel can be overwritten if a better measurement from galaxy proeprty exist.
    finearr["vel"] = finetree["vel"] #
    finearr["nout"]=nnza_all.a2b(finetree["nstep"],"nstep","nout")

    for field in ["id", "idx", "pos", "vel", "nstep", "nout"]:
        fields_interp.remove(field)


    lbt = tc.zred2gyr(nnza_all.a2b(finetree["nstep"],"nstep","zred"), z_now=0)
    lbt_cell = tc.zred2gyr(nnza_cell.a2b(mainarr["nstep"],"nstep","zred"), z_now=0)


    for mar in mainarr:
        finearr["pos"][finearr["nout"] == mar["nout"]] = mar["pos"]
        finearr["vel"][finearr["nout"] == mar["nout"]] = mar["vel"]


    for field in fields_interp:
        # Begining of merger
        if mainarr[field].ndim == 2:
            for i in range(3):
                if do_smooth:
                    r_p = smooth(mainarr[field][:,i],
                                 window_len=5,
                                 clip_tail_zeros=False)
                    finearr[field][:,i] = np.interp(lbt, lbt_cell, r_p)
                else:
                    r_p = mainarr[field][:,i]
                finearr[field][:,i] = np.interp(lbt, lbt_cell, mainarr[field][:,i])
        else:
            if do_smooth:
                r_p = smooth(mainarr[field],
                             window_len=5,
                             clip_tail_zeros=False) # odd number results in +1 element in the smoothed array.
            else:
                r_p = mainarr[field]
            finearr[field] = np.interp(lbt, lbt_cell, r_p)

    return finearr


def serialize(allresults, all_final_idxs,  nnza, nnza_cell,
              istep_max = 50,
              prg_dir="./",
              out_dir="./",
              nstep_too_short_main = 100):
    nouts = nnza_cell.nnza["nout"][:istep_max]
    print("Considering nouts: ", nouts)

    allresults = get_all_results(nouts, prg_dir=prg_dir, out_dir =result_dir)
    """
    For an unknown reason, some of galaxies are repeatedly analized, and found in the result_lambda pickle.
    Those repeatition must be removed from the begining.
    Until then, go through one more step to remove duplicates
    """
    all_sample_idxs=pickle.load(open(prg_dir + "all_sample_idxs.pickle", "rb"))

    serial_out_dir = out_base+"result_serial/"
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


    info = info.Info(nout=nouts[0])
    tc = cosmology.Timeconvert(info, zred_now=0)

    all_fid_ok=[]
    all_fidx_ok=[]


    for i, this_idx in enumerate(all_final_idxs):
        fname = prg_dir + "{}_adp.pickle".format(this_idx)
        if not os.path.isfile(fname):
            # dump_prgs    broken in the middle.
            continue
        #print(i, "IDX=",this_idx)
        adp = pickle.load(open(fname, "rb"))
        if min(adp[0][0]["nstep"][adp[0][0]["nstep"] > 0]) > nstep_too_short_main:
            print("Too short main tree. SKIP")
            continue
        # Append age to maintree and mainresult.
        lbt = tc.zred2gyr(nnza_all.a2b(adp[0][0]["nstep"],"nstep","zred"),z_now=0)
        adp[0][0] = recf.append_fields(adp[0][0], "time", lbt)
        max_step=len(allresults)
        this_gal = serialize_results.Serial_result(adp)
        cnt_merger=0
        bad_main=False
        for i, this_sats in enumerate(adp):
            nout_results=[]

            for sat in this_sats:
                sat_results=[]
                for ss in sat:
                    nout=nnza_all.step2out([ss["nstep"]]) # NOT nnza_cell.
                    #print(nout)
                    if nout in nouts:
                        #print(nout, ss["idx"])
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
                    this_gal.main_arr = serialize_results.galresult2rec(sat_results, is_main=True)
                    #print(len(this_gal.main_arr.nstep), this_gal.main_arr.nstep.ptp())
                    if len(this_gal.main_arr.nstep) <= this_gal.main_arr.nstep.ptp():
                        #bad_main=True
                        this_gal.main_arr = fill_main(this_gal.main_arr, nnza_cell)
                    this_gal.finearr = interpol_fine(this_gal, nnza_cell, nnza_all, do_smooth=True)
                    #print("BAD", bad_main)

                elif len(sat_results) > 0 and sat_results[0].mstar > 0.0:
                    #print("merger2")
                    this_gal.add_merger(sat_results, sat)
                    cnt_merger+=1
                    #this_gal.mergers.append(serialize_results.get_merger_props(this_gal.main_arr,
                    #                        serialize_results.galresult2rec(sat_results)))
            if bad_main:
                print("Bad main. Break")
                break
            #print(i)
            this_gal.data.append(nout_results)


        pickle.dump(this_gal, open(serial_out_dir+"serial_result{}.pickle".format(this_idx), "wb"))
        all_fidx_ok.append(this_idx)
        #break

    np.savetxt(serial_out_dir+"all_fidx_ok.txt", all_fidx_ok, fmt="%d")

    return [pickle.load(open(serial_out_dir+"serial_result{}.pickle".format(this_idx), "rb")) for this_idx in all_fidx_ok]
