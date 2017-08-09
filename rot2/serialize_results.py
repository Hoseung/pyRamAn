import numpy as np
import utils.match as mtc

def galresult2rec(sat):
    """
    Convert result of a tree (sat or main) into recarray.
    Having them in array is useful for calculating merger properties.
    """
    sat_data=np.recarray(len(sat), dtype=[("nstep", "<i4"),
                                          ("nout", "<i4"),
                                          ("id", "<i4"),
                                          ("idx", "<i4"),
                                          ("mstar", "<f8"),
                                          ("mgas", "<f8"),
                                          ("mgas_cold", "<f8"),
                                          ("reff", "<f8"),
                                          ("lgas", "<f8", (3)),
                                          ("pos", "<f8", (3)),
                                          ("vel", "<f8", (3)),
                                          ("lvec", "<f8", (3)),
                                          ("nvec", "<f8", (3))])
    sat_data.nout = np.array([ss.nout for ss in sat])
    sat_data.nstep = np.array([ss.nstep for ss in sat])
    sat_data.id = np.array([ss.id for ss in sat])
    sat_data.idx= np.array([ss.idx for ss in sat])
    sat_data.mstar= np.array([ss.mstar for ss in sat])
    sat_data.reff= np.array([ss.reff for ss in sat])
    sat_data.pos = np.array([(ss.xc, ss.yc, ss.zc) for ss in sat])
    sat_data.vel = np.array([(ss.vxc, ss.vyc, ss.vzc) for ss in sat])

    for i, ss in enumerate(sat):
        if hasattr(ss,"gas_results"):
            sat_data.mgas[i]= ss.gas_results["mgas_tot"]
            sat_data.mgas_cold[i]= ss.gas_results["mgas_cold"]
            sat_data.lgas[i]= ss.gas_results["Ln_gas"]
        if ss.lvec is not None:
            sat_data.lvec[i]= ss.lvec
            sat_data.nvec[i]= ss.nvec

    return sat_data


class Serial_result():
    """
    Only part of steps in AllDirectProgenitors(ADP) is calculated.
    So both ADP and 'results' must be stored.
    """
    def __init__(self, adp):
        self.maintree = adp[0][0]
        self.alltree = adp
        self.fidx = adp[0][0]["idx"][0]
        self.data = []
        self.mergers = []

    def add_merger(self, results, sattree):
        self.mergers.append(Merger(self.main_arr, results, sattree))

    def fill_bad_measure(self, sat):
        """
        Both main and satellites need this treatment.
        No reason this should bound to the Merger class.
        """
        # fill bad measurements of main
        # Note that I just throw away bad sat measurements.
        # Have a look at the old (paper1) ipynb.
        pass

    def get_main_counterpart(self, main_arr, rec_sat):
        return main_arr[mtc.match_list_ind(main_arr["nstep"], rec_sat["nstep"])]

    #def get_orbit_from_tree(self, merger):

    def get_merger_props(self, merger):
        """
            To be called after self.main_arr is assigned.
            merger.sat is a recarrays.

            Every gal has to have all the properties.
            None-good galaxies must have been filtered out beforehand.
        """
        merger.merger_arr = np.zeros(len(merger.sat), dtype=merger_props)
        main_part = self.get_main_counterpart(self.main_arr, merger.sat)

        rel_pos = main_part.pos - merger.sat.pos
        rel_vel = main_part.vel - merger.sat.vel
        Js=np.cross(rel_pos, rel_vel)
        j_orbital = Js[:]/np.sqrt(np.einsum('...i,...i', Js, Js))[:,None]

        # spin alignment
        merger.merger_arr["nstep"]=main_part.nstep
        merger.merger_arr["id_p"]=main_part.id
        merger.merger_arr["idx_p"]=main_part.idx
        merger.merger_arr["id_s"]=merger.sat.id
        merger.merger_arr["idx_s"]=merger.sat.idx
        merger.merger_arr["m_s"]=merger.sat.mstar
        merger.merger_arr["m_p"]=main_part.mstar
        merger.merger_arr["rel_pos"]=rel_pos
        merger.merger_arr["rel_vel"]=rel_vel
        merger.merger_arr["jorbit"]=j_orbital
        merger.merger_arr["dist"]=np.sqrt(np.einsum('...i,...i', rel_pos,rel_pos)) * 1e3 # in Kpc
        merger.merger_arr["reff_p"]=main_part.reff
        merger.merger_arr["reff_s"]=merger.sat.reff
        merger.merger_arr["spinmag"]=np.sqrt(np.einsum('...i,...i', merger.sat.lvec,merger.sat.lvec))
        merger.merger_arr["orbitang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_part.nvec, j_orbital))
        merger.merger_arr["spinang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_part.nvec, merger.sat.nvec))

    def cal_passages(self, p_min=0.5):
        """
        Identify orbits with minimum period of p_min.
        If the period is shorter than that, assume the end of merger is near,
        and everthing thereafter is the effect of merger / final coalescence.
        """


merger_props = [("nstep", '<i8'),
                ("nout", '<i8'),
                ("id_p", '<i8'),
                ("idx_p", '<i8'),
                ("id_s", '<i8'),
                ("idx_s", '<i8'),
                ("m_s", '<f8'),
                ("m_p", '<f8'),
                ("dist", "<f8"),
                ("reff_p","<f8"),
                ("reff_s","<f8"),
                ("spinang", '<f8'),
                ("spinmag", '<f8'),
                ("orbitang", '<f8'),
                ("rel_pos", '<f8', (3,)),
                ("rel_vel", '<f8', (3,)),
                ("jorbit", '<f8', (3,))]


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
    def __init__(self, rec_main, sat, sattree):
        self.main_idx=rec_main[0].idx
        # if main and sat are OK.
        #merger.count_mergers(sat)
        self.main = rec_main
        self.sattree = sattree
        self.sat = galresult2rec(sat)
        #merger.merger_arr = merger.get_merger_props(rec_main, merger.sat)
        #self.cal_passages()
