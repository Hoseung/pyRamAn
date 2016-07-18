# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:33:27 2015

1) convert HM/GM output
2) generate .cfg file
3) run consistant tree



@author: hoseung
"""

import os
#import numpy as np
import subprocess

from tree.tree_builder import convert_halo_list


def run(wdir ='./', nout_ini=None, is_gal=False,
        out_dir_g='GalaxyMaker/', out_dir_d='halo/'):
    base = os.path.abspath(wdir) + '/'
    print(base)
    cluster = base.split('/')[-2]
    if is_gal:
        if nout_ini is None:
            nout_ini=11
        out_dir = base + out_dir_g
    else:
        if nout_ini is None:
            nout_ini=7
        out_dir = base + out_dir_d


    convert_halo_list(nout_ini=nout_ini, nout_fi = 187,
                      base=base, out_dir=out_dir, is_gal=is_gal,
                      nmax_fracmax_ratio=1.01,
                      nmax_fracmax_ratio2=1.00,
                      frac_max_small=0.5)
        
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
    #is_gal = input("Galaxy? (y,n. default = n) \n")
    #if is_gal in ["YES", "Yes", "yes", "Y", "y", True]:
    #    is_gal = True
    #else:
    #    is_gal = False

    here = os.path.abspath('./')
    #nout_ini=int(input("not_ini=? (12)"))
    #if nout_ini == "":
    #    nout_ini = 12
    is_gal = False
    nout_ini = 20
    #run(wdir = here + "/" , is_gal=is_gal, nout_ini=nout_ini,
    #    out_dir_g="halo/")
    #clusters=["29176"]

    # 29828 -> 29830
    clusters=["01605", "04466", "05427", "10002", "14172",
              "17891", "24954", "28930", "29172", "35663",
              "36413", "36415", "49096", "06098", "07206", "39990"][-1]
def run_all(clusters, is_gal=True, nout_ini=20):
    for cluster in clusters:
        run(wdir = cluster + "/" , is_gal=is_gal, nout_ini=nout_ini)
