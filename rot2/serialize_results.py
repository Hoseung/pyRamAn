import numpy as np
import utils.match as mtc
from rot2.serialize_module import smooth

merger_props = [("nstep", '<i8'),
                ("nout", '<i8'),
                ("id_p", '<i8'),
                ("idx_p", '<i8'),
                ("id_s", '<i8'),
                ("idx_s", '<i8'),
                ("m_s", '<f8'),
                ("m_p", '<f8'),
                ("size_s", '<f8'),
                ("size_p", '<f8'),
                ("dist", "<f8"),
                ("reff_p","<f8"),
                ("reff_s","<f8"),
                ("spinang", '<f8'),
                ("spinmag", '<f8'),
                ("orbitang", '<f8'),
                ("rel_pos", '<f8', (3,)),
                ("rel_vel", '<f8', (3,)),
                ("jorbit", '<f8', (3,))]

gal_props = [ ("nstep", "<i4"),
              ("nout", "<i4"),
              ("zred", "<f8"),
              ("id", "<i4"),
              ("idx", "<i4"),
              ("mstar", "<f8"),
              ("mgas", "<f8"),
              ("mgas_cold", "<f8"),
              ("rgal", "<f8"),
              ("reff", "<f8"),
              ("lambda_r", "<f8"),
              ("vsig", "<f8"),
              ("sfr01", "<f8"),
              ("sfr05", "<f8"),
              ("sfr1", "<f8"),
              ("sfr_area", "<f8"),
              ("lgas", "<f8", (3)),
              ("pos", "<f8", (3)),
              ("vel", "<f8", (3)),
              ("lvec", "<f8", (3)),
              ("nvec", "<f8", (3))]

#fields_interp = ["reff", "rgal", "mstar", "mgas", "mgas_cold", "lgas", "lvec", "nvec"]


def galresult2rec(sat, is_main=False, fill_missing=True):
    """
    Convert result of a tree (sat or main) into recarray.
    Having them in array is useful for calculating merger properties.
    """
    if is_main:
        dtype = gal_props + [("P_tidal", "<f8")]
    else:
        dtype = gal_props

    main_arr=np.recarray(len(sat), dtype=dtype)
    main_arr.nout = np.array([ss.nout for ss in sat])
    main_arr.nstep = np.array([ss.nstep for ss in sat])
    main_arr.id = np.array([ss.id for ss in sat])
    main_arr.idx= np.array([ss.idx for ss in sat])
    main_arr.mstar= np.array([ss.mstar for ss in sat])
    main_arr.reff= np.array([ss.reff for ss in sat])
    main_arr.pos = np.array([(ss.xc, ss.yc, ss.zc) for ss in sat])
    main_arr.vel = np.array([(ss.vxc, ss.vyc, ss.vzc) for ss in sat])

    for i, ss in enumerate(sat):
        if hasattr(ss,"gas_results"):
            main_arr.mgas[i]= ss.gas_results["mgas_tot"]
            main_arr.mgas_cold[i]= ss.gas_results["mgas_cold"]
            main_arr.lgas[i]= ss.gas_results["Ln_gas"]
        if ss.lvec is not None:
            main_arr.lvec[i]= ss.lvec
            main_arr.nvec[i]= ss.nvec
            # If there is Lvec, there is rgal.
            main_arr.rgal[i]= ss.rgal
            try:
                main_arr.lambda_r[i]=ss.lambda_r[0]
            except:
                main_arr.lambda_r[i]=ss.lambda_r
            if is_main:
                if hasattr(ss, 'P'):
                    main_arr.P_tidal[i] = ss.P
                else:
                    main_arr.P_tidal[i] = 0
        try:
            main_arr["vsig"][i]=ss.vsig_results["V_sig"]
        except:
            main_arr["vsig"][i]=0

        if hasattr(ss, "sfr_results"):
            main_arr["sfr01"][i]=ss.sfr_results["sfrs"][0]
            main_arr["sfr05"][i]=ss.sfr_results["sfrs"][1]
            main_arr["sfr1"][i]=ss.sfr_results["sfrs"][2]
            main_arr["sfr_area"][i]=ss.sfr_results["area"]


    #if fill_missing and len(main_arr.nstep) <= main_arr.nstep.ptp():
    #    smooth_all(main_arr)

    return main_arr


def fill_main(mainarr, nnza_cell, tc):
    # Set up a new array.
    # [3:] because the lambda_reuslts are mixed-up.
    # remove later.
    new_nouts = nnza_cell.nnza["nout"][3:3+mainarr["nstep"].ptp()+1]
    newarr = np.zeros(len(new_nouts), dtype=mainarr.dtype)
    # It's easy to fill nouts and nsteps.
    newarr["nout"]=new_nouts
    newarr["nstep"]=nnza_cell.a2b(newarr["nout"], "nout", "nstep")

    interp_fields = list(mainarr.dtype.names)
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


def interpol_fine(this_gal, nnza_cell, nnza_all, tc, do_smooth=True):
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

    # Physical time based interpolation
    lbt = tc.zred2gyr(nnza_all.a2b(finetree["nstep"],"nstep","zred"), z_now=0)
    lbt_cell = tc.zred2gyr(nnza_cell.a2b(mainarr["nstep"],"nstep","zred"), z_now=0)

    for mar in mainarr:
        finearr["pos"][finearr["nout"] == mar["nout"]] = mar["pos"]
        finearr["vel"][finearr["nout"] == mar["nout"]] = mar["vel"]


    for field in fields_interp:
        # Begining of merger
        #print("Main tree interpol", field)
        #print(mainarr[field][15])
        if mainarr[field].ndim == 2:
            for i in range(3):
                if do_smooth:
                    r_p = smooth(mainarr[field][:,i],
                                 window_len=5,
                                 clip_tail_zeros=False)
                    #finearr[field][:,i] = np.interp(lbt, lbt_cell, r_p)
                else:
                    r_p = mainarr[field][:,i]
                #print("3", r_p)
                finearr[field][:,i] = np.interp(lbt, lbt_cell, r_p)
        else:
            if do_smooth:
                r_p = smooth(mainarr[field],
                             window_len=5,
                             clip_tail_zeros=False) # odd number results in +1 element in the smoothed array.
            else:
                r_p = mainarr[field]
            #print("1", r_p, lbt, lbt_cell)
            #dat =
            #print(dat)
            finearr[field] = np.interp(lbt, lbt_cell, r_p)

    return finearr



class Serial_result():
    """
    Only part of steps in AllDirectProgenitors(ADP) is calculated.
    So both ADP and 'results' must be stored.
    """
    def __init__(self, adp, nnza_all):
        self.set_maintree(adp[0][0], nnza_all)
        #self.maintree = adp[0][0]
        self.alltree = adp
        self.fidx = adp[0][0]["idx"][0]
        self.data = []
        self.mergers = []

    def set_maintree(self, tree, nnza_all):
        tree = tree[np.where(tree["nstep"] > 0)[0]]
        self.maintree = np.zeros(len(tree["nstep"]),
                                 dtype=gal_props + [("P_tidal", "<f8"),
                                                    ("time", "<f8")])
        self.maintree["pos"] = tree["xp"]
        self.maintree["vel"] = tree["vp"]
        self.maintree["id"] = tree["id"]#.compressed()
        self.maintree["idx"] = tree["idx"]#.compressed()
        self.maintree["nstep"] = tree["nstep"]#.compressed()
        self.maintree["mstar"] = tree["m"]#.compressed()
        self.maintree["zred"] = nnza_all.a2b(self.maintree["nstep"],"nstep","zred")
        self.maintree["nout"] = nnza_all.a2b(self.maintree["nstep"],"nstep","nout")
        #self.maintree["time"] = tree[""].compressed()

    def cal_fine_arr(self, lbt, lbt_cell, do_smooth=True):
        finetree=self.maintree
        mainarr = self.main_arr
        finearr = np.zeros(len(finetree),dtype=mainarr.dtype)

        #lbt = tc.zred2gyr(nnza_all.a2b(finetree["nstep"],"nstep","zred"), z_now=0)
        #lbt_cell = tc.zred2gyr(nnza_cell.a2b(mainarr["nstep"],"nstep","zred"), z_now=0)

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

        self.finearr = finearr

    def add_merger(self, results, sattree, nnza_all):
        #self.main_arr is assigned in serialize script.
        #print("min main_arr",self.main_arr["nout"])
        #sat_nouts = nnza_all.a2b(sattree["nstep"], "nstep", "nout")
        #print("All sattree", sat_nouts)

        # sat no longer than the main.
        #sattree = sattree[sattree["nstep"] > self.maintree["nstep"].min()]
        #mm = Merger(self.main_arr[self.main_arr["nout"] > min(sat_nouts)], results, sattree, self)
        # main longer than sat is not needed.
        # change field name.
        # nextsub -> nout
        name_org = sattree.dtype.names
        sattree.dtype.names = [dd if dd !="nextsub" else "nout" for dd in name_org]
        sattree["nout"] = nnza_all.a2b(sattree["nstep"], "nstep", "nout")
        mm = Merger(self.main_arr[self.main_arr["nout"] >= sattree["nout"].min()], results, sattree, self)
        mm.cal_merger_props()
        self.mergers.append(mm)

    def fill_bad_measure(self, sat):
        """
        Both main and satellites need this treatment.
        No reason this should bound to the Merger class.
        """
        # fill bad measurements of main
        # Note that I just throw away bad sat measurements.
        # Have a look at the old (paper1) ipynb.
        pass

    #def get_orbit_from_tree(self, merger):


    def cal_passages(self, p_min=0.5):
        """
        Identify orbits with minimum period of p_min.
        If the period is shorter than that, assume the end of merger is near,
        and everthing thereafter is the effect of merger / final coalescence.
        """

class Merger():
    """
        A class to hold satellite properties,
        and calculate merger related properties.

        1. Beginning and end of mergers
        2. Merger mass ratio
        3. total gas content
        4. Merger orbits
        5. Treatement of multiple mergers
        -> With Ultra-fine merger tree, the change of multiple mergers per snapshot is low.


        NOTE
        ----

    """
    def __init__(self, rec_main, sat, sattree, host_gal):
        self.main_idx=rec_main[0].idx
        # if main and sat are OK.
        #merger.count_mergers(sat)
        self.main = rec_main
        self.host = host_gal
        self.add_sattree(sattree)
        self.sat = galresult2rec(sat)

        self.cal_main_counterpart()
        #merger.merger_arr = merger.get_merger_props(rec_main, merger.sat)
        #self.cal_passages()

    def add_sattree(self, sattree):
        # sattree part earlier than the beginning of the main tree is not needed.
        #print("before", sattree["nstep"])
        #sattree=sattree[sattree["nstep"]>=min(self.host.maintree["nstep"][self.host.maintree["nstep"] > 0])]
        #print("After", sattree["nstep"])
        self.sattree = sattree

    def cal_main_counterpart(self):
        """
        Counterpart main result and main tree.
        """
        main_arr = self.host.main_arr
        # allow_swap = Fales
        # => Even if the sat is longer than the main,
        # indices for only common elements are returned.
        #print(main_arr["nstep"])
        imtc = mtc.match_list_ind(main_arr["nstep"], self.sat["nstep"],allow_swap=False)
        #print(len(imtc))
        self.main_part = main_arr[imtc]
        main_tree = self.host.maintree
        self.main_tree = main_tree[mtc.match_list_ind(main_tree["nstep"], self.sattree["nstep"],allow_swap=False)]
        #print("main",self.main_tree["nstep"], "sat",self.sat["nstep"])
        #return

    def cal_merger_props(self):
        """
            Calculate quantities often need for determining merger properties.
            Stores values in recarray and add to the given merger instance.

            Must be called after self.main_arr is assigned.

            merger.sat is a recarrays.

            Every gal has to have all the properties.
            None-good galaxies must have been filtered out beforehand.
        """
        self.merger_arr = np.zeros(len(self.sat), dtype=merger_props)
        main_part = self.main_part

        rel_pos = main_part["pos"] - self.sat.pos
        rel_vel = main_part["vel"] - self.sat.vel
        Js=np.cross(rel_pos, rel_vel)
        j_orbital = Js[:]/np.sqrt(np.einsum('...i,...i', Js, Js))[:,None]

        # spin alignment
        self.merger_arr["nstep"]=main_part["nstep"]
        self.merger_arr["id_p"]=main_part["id"]
        self.merger_arr["idx_p"]=main_part["idx"]
        self.merger_arr["id_s"]=self.sat["id"]
        self.merger_arr["idx_s"]=self.sat["idx"]
        self.merger_arr["m_s"]=self.sat["mstar"]
        self.merger_arr["m_p"]=main_part["mstar"]
        self.merger_arr["rel_pos"]=rel_pos
        self.merger_arr["rel_vel"]=rel_vel
        self.merger_arr["jorbit"]=j_orbital
        self.merger_arr["dist"]=np.sqrt(np.einsum('...i,...i', rel_pos,rel_pos)) * 1e3 # in Kpc
        self.merger_arr["reff_p"]=main_part["reff"]
        self.merger_arr["reff_s"]=self.sat.reff
        self.merger_arr["spinmag"]=np.sqrt(np.einsum('...i,...i', self.sat.lvec,self.sat.lvec))
        self.merger_arr["orbitang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_part["nvec"], j_orbital))
        self.merger_arr["spinang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_part["nvec"], self.sat.nvec))
