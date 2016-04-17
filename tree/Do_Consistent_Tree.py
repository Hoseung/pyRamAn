# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:33:27 2015

1) convert HM/GM output
2) generate .cfg file
3) run consistant tree



@author: hoseung
"""

import os
import numpy as np
import subprocess

def write_halo_bricks_ct(data, desc, nout, out_dir='./', is_gal=True):
    """
        format
        ------
        #ID DescID Mass Vmax Vrms Radius Rs Np X Y Z VX VY VZ JX JY JZ Spin
               
        ID: the halo ID, which must be unique across a single snapshot and must be at least 0.
        DescID: the halo's descendant id at the next timestep. If no descendant exists, this should be -1.
                At least 100 halos should have descendants (identified by particle-based methods) 
                in order for the consistency checks to work properly.
        Mass: the halo mass. This does not have to be Mvir, 
                but it must correspond to the mass with the radius in the Radius column. 
                The units for this must be Msun/h.
        Vmax: the maximum circular velocity, in units of km/s (physical, not comoving).
        Vrms: the velocity dispersion in units of km/s (physical, not comoving); 
                may be set to 0 if not available.
        Radius: the radius at which the mass is calculated in column 3. 
                The units must be in kpc / h. (Note that this is different from the position units!)
        Rs: the scale radius of the halo, in units of kpc / h; may be set to 0 if not available.
        Np: number of particles in the halo.
        X/Y/Z: 3D position of the halo, in units of Mpc / h.
                (Note that this is different from the radius units!)
        VX/VY/VZ: 3D velocity of the halo, in units of km/s
                (physical, not comoving).
        JX/JY/JZ: 3D angular momentum of the halo, in units of (Msun/h) * (Mpc/h) * km/s (physical, not comoving);
                may be set to 0 if not available.
        Spin: Dimensionless spin parameter; may be set to 0 if not available.
        
        
        output file name
        out_xyz.list (out_1.list, ..., out_100.list, ...)
        
        
        Note 1) 
            No Nans and Inf are allowed.
            
        Note 2) 
            You also need a list of snapshot numbers and aexps.
            
        Note 3)
            Works for both GalaxyMaker and HaloMaker.

    """
    def write_variables(f, variables):
        """
        Or, just data[['id', 'Orig_halo_id', ' ', and so on]] to access multiple fields at a time.
        
        """
        for variable in variables:
            f.write(str(variable) + " ")
        f.write("\n")
    

    if is_gal:
        field_mass = 'm'
        field_radius = 'r'
    else:
        field_mass = 'mvir'
        field_radius = 'rvir'            


    fn = out_dir + 'out_' + str(nout) + '.list'
    with open(fn, 'w') as f:
        f.write('#ID DescID Mass Vmax Vrms Radius Rs Np X Y Z VX VY VZ JX JY JZ Spin \n')
        for i in range(len(data['id'])):
            write_variables(f, (data['id'][i], desc[i], data[field_mass][i], \
            data['cvel'][i], 0 , data[field_radius][i], 0, data['np'][i], \
            data['x'][i], data['y'][i], data['z'][i],\
            data['vx'][i], data['vy'][i], data['vz'][i], \
            data['ax'][i], data['ay'][i], data['az'][i], \
            data['sp'][i]))
            
    f.close()
    
def write_scale(nouts, aexps, out_dir='./'):
    """
        fn = DescScales.txt
    """
    fn = out_dir + 'DescScales.txt'
    with open(fn, 'w') as f:
        for nout, aexp in zip(nouts, aexps):
            f.write(str(nout) + "  " + str(aexp) + "\n")
            
    f.close()


# convert HM/GM output
def convert_halo_list(nout_ini=11, nout_fi = 187, base='./', out_dir='./', is_gal=False):
    import load
    import tree.halomodule as hmo
    
    boxsize = 200
    nouts = range(nout_ini, nout_fi + 1)
    aexps = np.zeros(len(nouts))
    
    for inout, nout in enumerate(nouts):    
        if True:
            # If first snapshot, load data0.
            galcat0 = hmo.Halo(nout=nout, base=base,
                          halofinder="HM", is_gal=is_gal, return_id=True)
            galcat0.load()
            data0, idlists0 = galcat0.data, galcat0.idlists
        else:
            galcat0 = galcat1
            data0, idlists0 = galcat0.data, galcat0.idlists

        if nout == nouts[-1]:
            nhalo0 = galcat0.halnum + galcat0.subnum
            desc_list_final = np.zeros(nhalo0, dtype=int)
            desc_list_final.fill(-1)    
        else:       
            galcat1 = hmo.Halo(nout=nout+1, base=base,
                              halofinder="HM", is_gal=is_gal, return_id=True)
            galcat1.load()    
            data1, idlists1 = galcat1.data, galcat1.idlists
    
            # Some particles do not belong to any halo/galaxy. But most of them do. 
            # -> the length of hid_of_particls = number of total particles, not the sum of particles of each halo.
            # number of stellar particles change. 
            # -> Use larger value of the two snapshots. (always the next one?)
            npart = max([galcat0.nbodies,galcat1.nbodies]) + 1
            # + 1 so that particle id == array index.
    
            hid_of_particles0 = np.zeros(npart, dtype=int)
            hid_of_particles1 = np.zeros(npart, dtype=int)
            # idlist_per_halo to hid_of_particles
            for iha, ids in enumerate(idlists1):
                hid_of_particles1[ids] = data1['id'][iha]
    
            # progenitors of next halos
            # one halo can have multiple main progenitor.
            nhalo1 = galcat1.halnum + galcat1.subnum
            main_progenitors = [[]]* nhalo1
    
            # descendant of previous halos
            # one halo has only one descendants
            nhalo0 = galcat0.halnum + galcat0.subnum
            desc_list_final = np.zeros(nhalo0, dtype=int)
            for i, pid_per_gal in enumerate(idlists0):
                # Properties of current halos
                np_prog = data0['np'][i]
    
                # hid_of_particles1 == 0 means the particle belongs to no halo. 
                # Exclude them.
                particles_in_next_nout = np.bincount(hid_of_particles1[pid_per_gal])
                 #1 ~ 475
                all_desc_list = np.where(particles_in_next_nout > 0)[0][1:] # remove empty halos & zero-th halo
                if len(all_desc_list) == 0:
                    desc_list_final[i] = -1
                else:
                    # How many particles are coming from the progenitor?
                    np_received_list = particles_in_next_nout[all_desc_list]
                    i_max_np = np.argmax(np_received_list)
                    desc_np_max = all_desc_list[i_max_np]
                    # How much fraction of the descendant is coming from the progenitor?
                    #np_frac_list = [np / np_prog for np in np_received_list] # That is the fraction of prg.
                    frac_received_list = [np / np_desc for np, np_desc in zip(np_received_list, data1['np'][all_desc_list -1])]
    
                    i_max_frac = np.argmax(frac_received_list)
                    desc_frac_max = all_desc_list[i_max_frac]
    
                    if desc_np_max != desc_frac_max:
                        if np_received_list[i_max_np] > 10 * np_received_list[i_max_frac]:
                            # If one halo is significantly larger than the other, 
                            # which means much more fraction of the progenitor is inherited to the larger one,
                            # then the larger one is the descendant.
                            # By how much should it be larger?
    
                            desc_list_final[i] = frac_received_list[i_max_np]
                            # and, the other halo may or may not find it's progenitor.
                            # If the fraction if very high, it is very likely that the halo is coming out of nothing.
                            # Or at least the "real progenitor" had been misstaken as a part of another larger halo in the previous snapshot.
                        elif frac_received_list[i_max_frac] < 0.3:
                            if np_received_list[i_max_np] > 2 * np_received_list[i_max_frac]:
                                desc_list_final[i] = frac_received_list[i_max_np]
                            else:
                                desc_list_final[i] = frac_received_list[i_max_frac]
                    else:
                        desc_list_final[i] = desc_np_max
                        #print(i, desc_np_max)
    
        
        # Convert variables in appropriate units
        # data['cve'] is the circular velocit at Rvir. This is not the maximum circular velocity! 
        # angular momentum in 10**11 Msun * km/s * Mpc
        info = load.info.Info(nout=nout, base=base, load=True)
        data0['rvir'] *= boxsize * 1000 # in kpc/h
        data0['r'] *= boxsize * 1000 # in kpc/h
        data0['mvir'] = data0['mvir'] / info.h * 100 # Msun/h
        data0['m'] = data0['m'] / info.h * 100 # Msun/h
        data0['x'] *= boxsize
        data0['y'] *= boxsize
        data0['z'] *= boxsize
        data0['ax'] = data0['ax'] * 1e11/ info.h * 100 * 1000 / info.h * 100
        data0['ay'] = data0['ay'] * 1e11/ info.h * 100 * 1000 / info.h * 100
        data0['az'] = data0['az'] * 1e11/ info.h * 100 * 1000 / info.h * 100
        write_halo_bricks_ct(data0, desc_list_final, nout, out_dir = out_dir, is_gal=is_gal)
        aexps[inout] = info.aexp
    
    print("Write Descscale")    
    write_scale(nouts, aexps, out_dir = out_dir)



def run(wdir ='./', nout_ini=None, is_gal=False):
    base = os.path.abspath(wdir) + '/'
    cluster = base.split('/')[-2]
    if is_gal:
        if nout_ini is None:
            nout_ini=11
        out_dir = base + 'GalaxyMaker/'   
    else:
        if nout_ini is None:
            nout_ini=7
        out_dir = base + 'halo/'
    convert_halo_list(nout_ini=nout_ini, nout_fi = 187, base=base, out_dir=out_dir, is_gal=is_gal)
        
    #configure file template in consistent tree directory. or repo.
    f_cfg_template = '/home/hopung/Work/pyclusterevol/repo/CTree_template.cfg'
    ct_install_dir = "/home/hopung/Work/consistent_trees-1.01"
    
    
    
    for path in ['Outputs', 'Trees']:
        if not os.path.isdir(out_dir + path):
            os.mkdir(out_dir + path)
    
    
    # ConsistentTrees deals with multi-snaphsot linking taking all the halo properties into consideration.
    # So just find the descendants.
    # It is suffice to determine 'descendants' in a simple, classic way:
    # The halo that gets the most particles from the previous halo is the 'main' descendant. 
    # One possible improvement is to consider mass fraction, which is useful in case of halos
    # undergoing tidal stripping in cluster environments.
    # i.e., when a halo loses more than half of its particles to its host.
    
    # Sometimes.... 
    # 
    # Not matched
    # [17517   430   497]
    # [0.92151086327529064, 0.28067885117493474, 1.0]
    # 
    # one halo got most of the particles from the progenitor,
    # while another halo got 100% of particles for the progenitor.
    # But that 497-particled halo might be a temporary thing. 
    
    
    #%%
    # generate .cfg file
    inbase = out_dir
    scalefile = inbase + "DescScales.txt"
    outbase = inbase + 'Outputs'
    tree_outbase = inbase + 'Trees'
    hlist_outbase = inbase + 'Trees'
    
    if is_gal:
        mass_res_ok = '1e10'
        fname_cfg = inbase + cluster + '_GM.cfg'
    elif not is_gal: 
        mass_res_ok = '5e10'
        fname_cfg = inbase + cluster + '_HM.cfg'
    
    
    modify_strs = ["SCALEFILE", "INBASE", "OUTBASE", "TREE_OUTBASE", "HLIST_OUTBASE"]
    modify_vals=["MASS_RES_OK"]
    new_strs = [scalefile, inbase, outbase, tree_outbase, hlist_outbase]
    new_vals = [mass_res_ok]
    
    
    f_temp = open(f_cfg_template, 'r')     
    f_write = open(fname_cfg, 'w')
    
    for line in f_temp.readlines():
        start_with = line.split("=")[0].strip()
        if start_with in modify_strs:
            ind = modify_strs.index(start_with.strip())
            f_write.write(start_with + '= "' + new_strs[ind] + '"\n')
        elif start_with in modify_vals:
            ind = modify_vals.index(start_with.strip())
            f_write.write(start_with + '= ' + new_vals[ind] + '\n')
        else:
            f_write.write(line)
                
    f_write.close()
    f_temp.close()     
    #%%
    # run Consistent Trees
    os.chdir(ct_install_dir)
    subprocess.call(["perl", "do_merger_tree_np.pl", fname_cfg])
    subprocess.call(["perl", "halo_trees_to_catalog.pl", fname_cfg])
     

#%%   
if __name__ == "__main__":
#    clusters=["01605", "04466", "05427", "074010", "10002", "14172",
#              "17891", "24954", "28928", "29172", "29176", "35663",
#              "36413", "36415_000", "49096", "6098", "7206", "39990"]
    
    is_gal = input("Galaxy? (y,n. default = n) \n")
    if is_gal in ["YES", "Yes", "yes", "Y", "y", True]:
        is_gal = True
    else:
        is_gal = False

    here = os.path.abspath('./')
#    for cluster in clusters[6:7]:
#    print("Processing....")
#     print(cluster)
    nout_ini=int(input("not_ini=? (12)"))
    if nout_ini == "":
        nout_ini = 12
    run(wdir = here + "/" , is_gal=is_gal, nout_ini=nout_ini)
