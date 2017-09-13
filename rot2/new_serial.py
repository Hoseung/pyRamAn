from utils import match as mtc
import numpy as np

dtype_maintree = [('nstep', '<i4'), ('nout', '<i4'), ('lbt', '<f8'),
                  ('idx', '<i4'), ('id', '<i4'),
                  ('m', '<f8'),
                  ('pos', '<f8', (3,)), ('vel', '<f8', (3,)), ('rgal', '<f8'),
                  ('cvel', '<f8'),
                  ('P_tidal', '<f8'),
                  ('d5', '<f8'),
                  ('d10', '<f8'),
                  ('d50', '<f8')]

dtype_merger_props = [("nstep", '<i8'),
                ("nout", '<i8'),
                ("id", '<i8'),
                ("idx", '<i8'),
                ("dist", "<f8"),
                ("orbitang", '<f8'),
                ("m", "<f8"),
                ("rel_pos", '<f8', (3,)),
                ("rel_vel", '<f8', (3,)),
                ("jorbit", '<f8', (3,)),
                # Below are to be interpolated from measurements.
                ("spinang", '<f8'),
                ("spinmag", '<f8'),
                ("mstar", '<f8'),
                ("m_frac", '<f8'),
                ("m_gas", '<f8'),
                ("m_gas_cold", '<f8'),
                ("size_s", '<f8'),
                ("size_p", '<f8'),
                ("rgal_sum","<f8"),
                ("reff_s","<f8")]

dtype_main_results = [("nstep", "<i4"),
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
                      ("eps", "<f8"),
                      ("vsig", "<f8"),
                      ("sfr01", "<f8"),
                      ("sfr05", "<f8"),
                      ("sfr1", "<f8"),
                      ("sfr_area", "<f8"),
                      ("lgas", "<f8", (3)),
                      ("pos", "<f8", (3)),
                      ("vel", "<f8", (3)),
                      ("lvec", "<f8", (3)),
                      ("nvec", "<f8", (3)),
                      ("P_tidal_h", "<f8")]

dtype_sat_results = [ ("nout", "<i4"),
                      ("zred", "<f8"),
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


class Maingal():
    def __init__(self, adp, nnza_all, nnouts, nouts_cell):
        self.i_last_main=0 # main_arr data input index
        self.fidx = adp[0][0]["idx"][0]
        # maintree is poppeds
        self.set_maintree(adp[0].pop(0), nnza_all)
        # Now Only sats are remaing.
        self.main_data=np.zeros(nnouts, dtype=dtype_main_results)
        self.main_data["nout"] = nouts_cell[:nnouts]
        self.mergers=[]
        self.sat_data=[]
        self.add_mergers_tree(adp, nouts_cell)
        self.n_minor = 0
        self.n_major = 0

    def set_maintree(self, tree, nnza_all):
        tree = tree[np.where(tree["nstep"] > 0)[0]]
        self.maintree = np.zeros(len(tree["nstep"]), dtype=dtype_maintree)
        self.maintree["pos"] = tree["xp"]
        self.maintree["vel"] = tree["vp"]
        self.maintree["id"] = tree["id"]
        self.maintree["idx"] = tree["idx"]
        self.maintree["nstep"] = tree["nstep"]
        self.maintree["m"] = tree["m"] # to be
        self.maintree["lbt"] = nnza_all.a2b(self.maintree["nstep"],"nstep","lbt")
        self.maintree["nout"] = nnza_all.a2b(self.maintree["nstep"],"nstep","nout")

    def add_main_result(self, ss):
        # 63 steps
        il=self.i_last_main
        self.measurements["nout"][il] = ss.nout
        self.measurements["nstep"][il] = ss.nstep
        self.measurements["id"][il] = ss.id
        self.measurements["idx"][il] = ss.idx
        self.measurements["mstar"][il] = ss.mstar
        self.measurements["reff"][il] = ss.reff
        self.measurements["pos"][il,0] = ss.xc
        self.measurements["pos"][il,1] = ss.yc
        self.measurements["pos"][il,2] = ss.zc
        self.measurements["vel"][il,0] = ss.vxc
        self.measurements["vel"][il,1] = ss.vyc
        self.measurements["vel"][il,2] = ss.vzc

        self.measurements["lambda_r"][il] = data.nout

        if hasattr(ss, "gas_results"):
            self.measurements["mgas"] = ss.gas_results["mgas_tot"]
            self.measurements["mgas_col"] = ss.gas_results["mgas_cold"]
            self.measurements["lgas"]= ss.gas_results["Ln_gas"]
        if ss.lvec is not None:
            self.measurements["lvec"][il] = ss.lvec
            self.measurements["nvec"][il] = ss.nvec
            # If there is Lvec, there is rgal.
            self.measurements["rgal"][il] = ss.rgal
            try:
                self.measurements["lambda_r"][il]=ss.lambda_r[0]
            except:
                self.measurements["lambda_r"][il]=ss.lambda_r

        try:
            self.measurements["vsig"][il]=ss.vsig_results["V_sig"]
        except:
            pass

        if hasattr(ss, "sfr_results"):
            self.measurements["sfr01"][il]=ss.sfr_results["sfrs"][0]
            self.measurements["sfr05"][il]=ss.sfr_results["sfrs"][1]
            self.measurements["sfr1"][il]=ss.sfr_results["sfrs"][2]
            self.measurements["sfr_area"][il]=ss.sfr_results["area"]
    #def finalize_main_result

    def add_mergers_tree(self, adp, nouts_cell):
        """
            Add merger_arr from sat fine merger trees.
        """
        # sat should be no longer than the main.
        #sattree = sattree[sattree["nstep"] > self.maintree["nstep"].min()]
        #mm = Merger(self.main_arr[self.main_arr["nout"] > min(sat_nouts)], results, sattree, self)
        # main longer than sat is not needed.
        # change field name.
        # nextsub -> nout
        for adps_now in adp:
            for sattree in adps_now:
                marr = cal_merger_props_tree(self.maintree, sattree)
                if marr is not None:
                    self.mergers.append(marr)

        # sat_data to hold sat measurements
        for mm in self.mergers:
            matched_nouts =np.intersect1d(mm["nout"],nouts_cell)
            sd = np.zeros(len(matched_nouts), dtype=dtype_sat_results)
            sd["nout"]=matched_nouts[::-1]
            self.sat_data.append(sd)


        #name_org = sattree.dtype.names
        #sattree.dtype.names = [dd if dd !="nextsub" else "nout" for dd in name_org]
        #sattree["nout"] = nnza_all.a2b(sattree["nstep"], "nstep", "nout")
        #mm = Merger(self.main_arr[self.main_arr["nout"] >= sattree["nout"].min()], results, sattree, self)

def cal_merger_props_tree(maintree, sattree,
                           mass_ratio_cut=1/50, min_len=15):
    """
        Calculate merger related quantities, this time with only tree information.

        Parameters
        ----------
        mass_ratio_cut: 1/50
            If the maximum mass ratio is below mass_ratio_cut, ignore.
        min_len : 15
            Sat tree shorter than 15 is ignored.
            (15 in fine tree ~= 1 in coarse tree)

        NOTE
        ----
            More through mass ratio test must be done onnce
            the beginning of merger is determined.
    """
    if len(sattree) < 15:
        return None

    # only overlapping part of trees.
    common_nsteps = np.intersect1d(maintree["nstep"], sattree["nstep"])[::-1]
    if len(common_nsteps) < 15:
        return None

    main_tree_part = maintree[mtc.match_list_ind(maintree["nstep"],common_nsteps)]
    sattree = sattree[mtc.match_list_ind(sattree["nstep"],common_nsteps)]

    assert (np.all(main_tree_part["nstep"] == sattree["nstep"])), "main and tree are incompatible"

    # If the sat is too small, throw away.
    if np.max(sattree["m"]/main_tree_part["m"]) < mass_ratio_cut:
        return None

    rel_pos = main_tree_part["pos"] - sattree["xp"]
    # if more than half is confused with the main tree,
    # throw it away.
    if np.sum(rel_pos[:,0] == 0) > len(rel_pos)/2:
        return None
    rel_vel = main_tree_part["vel"] - sattree["vp"]

    merger_arr = np.zeros(len(sattree), dtype=dtype_merger_props)
    merger_arr["jorbit"]=np.cross(rel_pos, rel_vel)

    # spin alignment
    merger_arr["nstep"]=sattree["nstep"]
    merger_arr["nout"]=main_tree_part["nout"]
    merger_arr["id"]=sattree["id"]
    merger_arr["m"] = sattree["m"] # Just for checking
    merger_arr["idx"]=sattree["idx"]
    merger_arr["rel_pos"]=rel_pos
    merger_arr["rel_vel"]=rel_vel

    merger_arr["dist"]=np.sqrt(np.einsum('...i,...i', rel_pos,rel_pos)) * 1e3 # in Kpc
    return merger_arr

def cal_merger_props_measurments(main_arr, sat_arr):
    merger_arr["orbitang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_tree_part["nvec"], j_orbital))
    merger_arr["reff_p"]=main_tree_part["reff"]
    merger_arr["reff_s"]=self.sat.reff
    merger_arr["spinmag"]=np.sqrt(np.einsum('...i,...i', self.sat.lvec,self.sat.lvec))
    merger_arr["spinang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_tree_part["nvec"], self.sat.nvec))
    return merger_arr
