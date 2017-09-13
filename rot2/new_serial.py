from utils import match as mtc

dtype_maintree = [('nstep', '<i4'), ('nout', '<i4'), ('lbt', '<f8'),
                  ('idx', '<i4'), ('id', '<i4'),
                  ('mstar', '<f8'),
                  ('pos', '<f8', (3,)), ('vel', '<f8', (3,)), ('rgal', '<f8'),
                  ('cvel', '<f8'),
                  ('P_tidal', '<f8'),
                  ('d5', '<f8'),
                  ('d10', '<f8'),
                  ('d50', '<f8')]

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



class Maingal():
    def __init__(self, adp, nnza_all):
        self.i_last_main=0 # main_arr data input index
        self.fidx = adp[0][0]["idx"][0]
        # maintree is poppeds
        self.set_maintree(adp[0].pop(0), nnza_all)
        # Now Only sats are remaing.
        self.measurements=np.zeros(nnouts, dtype=dtype_main_results)
        self.sats=[]
        self.add_sats(adp)

    def set_maintree(self, tree, nnza_all):
        tree = tree[np.where(tree["nstep"] > 0)[0]]
        self.maintree = np.zeros(len(tree["nstep"]), dtype=dtype_maintree)
        self.maintree["pos"] = tree["xp"]
        self.maintree["vel"] = tree["vp"]
        self.maintree["id"] = tree["id"]
        self.maintree["idx"] = tree["idx"]
        self.maintree["nstep"] = tree["nstep"]
        self.maintree["mstar"] = tree["m"]
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

    def add_sats(self, adp):

    def add_merger(self, results, sattree, nnza_all):
        # sat should be no longer than the main.
        #sattree = sattree[sattree["nstep"] > self.maintree["nstep"].min()]
        #mm = Merger(self.main_arr[self.main_arr["nout"] > min(sat_nouts)], results, sattree, self)
        # main longer than sat is not needed.
        # change field name.
        # nextsub -> nout
        name_org = sattree.dtype.names
        sattree.dtype.names = [dd if dd !="nextsub" else "nout" for dd in name_org]
        sattree["nout"] = nnza_all.a2b(sattree["nstep"], "nstep", "nout")
        mm = Merger(self.main_arr[self.main_arr["nout"] >= sattree["nout"].min()], results, sattree, self)

        self.mergers.append(cal_merger_props(self.maintree, sattree))
        self.mergers.append(mm)


def cal_merger_props_trees(maintree, sattree):
    """
        Calculate quantities often need for determining merger properties.
        Stores values in recarray and add to the given merger instance.

        Must be called after self.main_arr is assigned.

        merger.sat is a recarrays.

        Every gal has to have all the properties.
        No good galaxies must have been filtered out beforehand.
    """

    main_part = maintree[mtc.match_list_ind(maintree["nstep"],sattree["nstep"])]
    assert (np.all(main_part["nstep"] == sattree["nstep"])), "main and tree are incompatible"

    merger_arr = np.zeros(len(self.sat), dtype=merger_props)

    rel_pos =main_part["pos"] - sattree.pos
    rel_vel = main_part["vel"] - sattree.vel
    Js=np.cross(rel_pos, rel_vel)
    j_orbital = Js[:]/np.sqrt(np.einsum('...i,...i', Js, Js))[:,None]

    # spin alignment
    merger_arr["nstep"]=sattree["nstep"]
    merger_arr["id"]=sattree["id"]
    merger_arr["idx"]=sattree["idx"]
    merger_arr["rel_pos"]=rel_pos
    merger_arr["rel_vel"]=rel_vel

    merger_arr["dist"]=np.sqrt(np.einsum('...i,...i', rel_pos,rel_pos)) * 1e3 # in Kpc
    merger_arr["orbitang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_part["nvec"], j_orbital))
    return merger_arr

def cal_merger_props_measurments(main_arr, sat_arr):
    merger_arr["jorbit"]=j_orbital
    merger_arr["reff_p"]=main_part["reff"]
    merger_arr["reff_s"]=self.sat.reff
    merger_arr["spinmag"]=np.sqrt(np.einsum('...i,...i', self.sat.lvec,self.sat.lvec))
    merger_arr["spinang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_part["nvec"], self.sat.nvec))
    return merger_arr
