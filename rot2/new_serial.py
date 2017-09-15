from utils import match as mtc
import numpy as np

dtype_finedata = [('nstep', '<i4'),
                  ('nout', '<i4'),
                  ('lbt', '<f8'),
                  ('idx', '<i4'),
                  ('id', '<i4'),
                  ('mstar', '<f8'),
                  ('pos', '<f8', (3,)),
                  ('vel', '<f8', (3,)),
                  ('rgal', '<f8'),
                  ('cvel', '<f8'),
                  ('P_tidal', '<f8'),
                  ("P_tidal_h", "<f8"),
                  ('d5', '<f8'),
                  ('d10', '<f8'),
                  ('d50', '<f8'),
                  ("mgas", "<f8"),
                  ("mgas_cold", "<f8"),
                  ("reff", "<f8"),
                  ("lambda_r", "<f8"),
                  ("eps", "<f8"),
                  ("vsig", "<f8"),
                  ("sfr01", "<f8"),
                  ("sfr05", "<f8"),
                  ("sfr1", "<f8"),
                  ("sfr_area", "<f8"),
                  ("lgas", "<f8", (3)),
                  ("lvec", "<f8", (3)),
                  ("nvec", "<f8", (3))]

dtype_merger_props = [("lbt", '<f8'),
                    ("nout", '<i8'),
                    ("id", '<i8'),
                    ("idx", '<i8'),
                    ("dist", "<f8"),
                    ("orbitang", '<f8'),
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
                    ("rgal", '<f8'),
                    ("rgal_sum","<f8"),
                    ("reff_s","<f8")]

dtype_main_results = [("nstep", "<i4"),
                      ("nout", "<i4"),
                      ("lbt", "<f8"),
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
                      ("lbt", "<f8"),
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

min_sattree_length = 21

class Maingal():
    def __init__(self, adp, nnza_all, nnouts, nnza_cell):
        self.i_last_main=0 # main_arr data input index
        self.fidx = adp[0][0]["idx"][0]
        # finedata is poppeds
        self.set_finedata(adp[0].pop(0), nnza_all)

        # Now Only sats are remaing.
        self.main_data=np.zeros(nnouts, dtype=dtype_main_results)
        self.main_data["nout"] = nnza_cell.nnza["nout"][:nnouts]
        self.main_data["lbt"] = nnza_cell.nnza["lbt"][:nnouts]
        # What if main gal measurements are not available at early nouts?
        self.mergers=[]
        self.sat_data=[]
        self.add_mergers_tree(adp, nnza_cell, nnouts)
        self.n_minor = 0
        self.n_major = 0
        self.nout_data = nnza_cell.nnza["nout"][:nnouts]

    def set_finedata(self, tree, nnza_all):
        tree = tree[np.where(tree["nstep"] > 0)[0]]
        self.finedata = np.zeros(len(tree["nstep"]), dtype=dtype_finedata)
        self.finedata["pos"] = tree["xp"]
        self.finedata["vel"] = tree["vp"]
        self.finedata["id"] = tree["id"]
        self.finedata["idx"] = tree["idx"]
        self.finedata["nstep"] = tree["nstep"]
        self.finedata["mstar"] = tree["m"] # to be
        self.finedata["lbt"] = nnza_all.a2b(self.finedata["nstep"],"nstep","lbt")
        self.finedata["nout"] = nnza_all.a2b(self.finedata["nstep"],"nstep","nout")


    def add_mergers_tree(self, adp, nnza_cell, nnouts):
        """
            Add merger_arr from sat fine merger trees.
        """
        # sat should be no longer than the main.
        #sattree = sattree[sattree["nstep"] > self.finedata["nstep"].min()]
        #mm = Merger(self.main_arr[self.main_arr["nout"] > min(sat_nouts)], results, sattree, self)
        # main longer than sat is not needed.
        # change field name.
        # nextsub -> nout
        for adps_now in adp:
            for sattree in adps_now:
                marr = cal_merger_props_tree(self, sattree)
                if marr is not None:
                    self.mergers.append(marr)

        # sat_data to hold sat measurements
        for mm in self.mergers:
            matched_nouts =np.intersect1d(mm["nout"],nnza_cell.nnza["nout"][:nnouts])
            sd = np.full(len(matched_nouts), np.nan, dtype=dtype_sat_results)
            sd["nout"]=matched_nouts[::-1]
            sd["lbt"]=nnza_cell.a2b(sd["nout"], "nout", "lbt")
            self.sat_data.append(sd)


def cal_merger_props_tree(this_gal, sattree,
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
    finedata = this_gal.finedata
    common_nsteps = np.intersect1d(finedata["nstep"], sattree["nstep"])[::-1]
    if len(common_nsteps) < min_sattree_length:
        return None

    main_tree_part = finedata[mtc.match_list_ind(finedata["nstep"],common_nsteps)]
    sattree = sattree[mtc.match_list_ind(sattree["nstep"],common_nsteps)]
    assert (np.all(main_tree_part["nstep"] == sattree["nstep"])), "main and tree are incompatible"

    # Even the trees are long and overlapping, measured part can only minimally overlap.
    # Abort in such cases.
    if this_gal.main_data["nout"][this_gal.main_data["nout"] > 0][-3] > main_tree_part["nout"].max():
        return None

    # If the sat is too small, throw away.
    if np.max(sattree["m"]/main_tree_part["mstar"]) < mass_ratio_cut:
        return None

    rel_pos = main_tree_part["pos"] - sattree["xp"]
    # if more than half is confused with the main tree,
    # throw it away.
    if np.sum(rel_pos[:,0] == 0) > len(rel_pos)/2:
        return None

    rel_vel = main_tree_part["vel"] - sattree["vp"]

    merger_arr = np.zeros(len(sattree), dtype=dtype_merger_props)
    merger_arr["nout"]=main_tree_part["nout"]
    # spin alignment
    #merger_arr["nstep"]=sattree["nstep"]

    merger_arr["id"]=sattree["id"]
    merger_arr["mstar"] = sattree["m"] # Just for checking
    merger_arr["idx"]=sattree["idx"]
    merger_arr["rel_pos"]=rel_pos
    merger_arr["rel_vel"]=rel_vel
    return merger_arr

def cal_merger_props_measurments(main_arr, sat_arr):
    merger_arr["orbitang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_tree_part["nvec"], j_orbital))
    merger_arr["reff_p"]=main_tree_part["reff"]
    merger_arr["reff_s"]=self.sat.reff
    merger_arr["spinmag"]=np.sqrt(np.einsum('...i,...i', self.sat.lvec,self.sat.lvec))
    merger_arr["spinang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_tree_part["nvec"], self.sat.nvec))
    return merger_arr



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
        # Todo
        # suppress nvec being nan when measuring it from stars.
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
