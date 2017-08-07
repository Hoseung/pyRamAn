import numpy as np
import utils.match as mtc

class Serial_result():
    def __init__(self, adp):
        self.maintree = adp[0][0]
        #self.prg_ids = aprg[1]
        self.alltree = adp
        self.fidx = adp[0][0]["idx"][0]
        self.mergers = []
        self.data = []
        self.sat_data = None


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



def galresult2rec(sat):
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



def get_merger_props(rec_main, rec_sat):
    """
        main, sat are recarrays.

        Every gal has to have all the properties.
        None-good galaxies must have been filtered out beforehand.
    """
    merger_list = np.zeros(len(rec_sat), dtype=merger_props)
    main_part = rec_main[mtc.match_list_ind(rec_main.nstep, rec_sat.nstep)]

    rel_pos = main_part.pos - rec_sat.pos
    rel_vel = main_part.vel - rec_sat.vel
    Js=np.cross(rel_pos, rel_vel)
    j_orbital = Js[:]/np.sqrt(np.einsum('...i,...i', Js, Js))[:,None]

    # spin alignment
    merger_list["nstep"]=main_part.nstep
    #merger_list["nout"]=main_part.nout
    merger_list["id_p"]=main_part.id
    merger_list["idx_p"]=main_part.idx
    merger_list["id_s"]=rec_sat.id
    merger_list["idx_s"]=rec_sat.idx
    merger_list["m_s"]=rec_sat.mstar
    merger_list["m_p"]=main_part.mstar
    merger_list["rel_pos"]=rel_pos
    merger_list["rel_vel"]=rel_vel
    merger_list["jorbit"]=j_orbital
    merger_list["dist"]=np.sqrt(np.einsum('...i,...i', rel_pos,rel_pos)) * 1e3 # in Kpc
    merger_list["reff_p"]=main_part.reff
    merger_list["reff_s"]=rec_sat.reff
    merger_list["spinmag"]=np.sqrt(np.einsum('...i,...i', rec_sat.lvec,rec_sat.lvec))
    merger_list["orbitang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_part.nvec, j_orbital))
    merger_list["spinang"] = 180./np.pi*np.arccos(np.einsum('...i,...i',main_part.nvec, rec_sat.nvec))

    return merger_list

class Mergers():
    """
        A class to hold satellite properties of a main galaxy tree,
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
    def __init__(self, main, sat, sattree):
        self.main_idx=main[0].idx
        # if main and sat are OK.
        self.count_mergers(sat)
        self.main = main
        self.sat = sat
        self.get_mergers(main, sat)

    def count_mergers(self, sat):
        cnt = 0
        for ss in sat:
            cnt += len(ss)
        self.cnt=cnt

    def fill_bad_measure_main(self, main):
        """
        Both main and satellites need this treatment.
        No reason this should bound to the Mergers class. 
        """
        # fill bad measurements of main
        # Note that I just throw away bad sat measurements.
        # Have a look at the old (paper1) ipynb.
        pass

    def get_mergers(self, main, sat):
        self.mergerlist = np.zeros(self.cnt, dtype=merger_props)
        i_arr=0
        for i, (this_main, this_sats) in enumerate(zip(main, sat)):
            if this_sats is not None:
                if len(this_sats) > 0:
                    # There are mergers.
                    for this_sat in this_sats:
                        if this_sat == -1:
                            self.mergerlist["sid"]=-1
                        elif hasattr(this_main,"nout"):
                            # The main galaxy may have poor measurement at a certain snapshot.
                            # skip and fill the hole by smoothng later.
                            get_merger_props(this_main, this_sat, self.mergerlist[i_arr])
                        i_arr+=1
