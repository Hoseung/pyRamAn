# -*- coding: utf-8 -*-
"""
Created on Fri May 29 02:43:32 2015

save copy of rockstar output files in ./new directory with 
corrected snapshot counts.

For example, if you ran rockstar for 090 - 099 snapshots:

rockstar_halos/halos_0.1.ascii -> rockstar_halos/halos_90.1.ascii
rockstar_halos/out_0.list -> rockstar_halos/out_90.list


It refers to aout list in .nml file to guess correct nout value. 

Pickle Rockstar ASCII outputs. ??

@author: hoseung
"""

from os.path import isdir
from os import path
from os import mkdir
from shutil import copyfile
import pickle
from glob import glob


def aexp2nout(nml, aexp):
    """ Accepts a namelist instance and an aexp(float).
        Returns closest aexp from the nml, and nout for that aexp
        
        Notes
        -----
        0.85 is represented as 0.849999999997. is that fine?
        >>> np.round(0.849999999997, decimals=2)
        >>> 0.84999999999999998
        >>> np.round(0.849999999997, decimals=2) == 0.85
        >>> True
        
        but, 
        >>> nml.aexp[155] == 0.85
        >>> False
        
        Hmmm....
    """
    import numpy as np
    d = np.abs(nml.aout - aexp)
    aout = nml.aout[d.argmin()]
    
    return aout, d.argmin() + 1

def run(base = './',
        rs_dir = 'rockstar_halos2/',
        dump_pickle=True, 
        nout_list=None,
        fnml=None):
    
    # should not include "halos_1" string.
    # this is used to search for output files. (rockstar_halos1/halos_1.0.ascii)
    
    new_dir = base + 'new/'   
    rsbase = base + rs_dir 
    # directory that pickled halos will go
    dir_out = path.normpath(base + 'halos_py')
    if not isdir(dir_out):
        mkdir(dir_out)
    
    print("go")
    
    # directory that re-named files go
    if not isdir(new_dir):
        mkdir(new_dir)
    
    if nout_list is None:
        houtput_list = glob(rsbase + 'halos_*.0.ascii')
    else:
        houtput_list = []
        for nout in nout_list:
            houtput_list.append(glob(rsbase + 'halos_' + str(nout) +'.0.ascii')[0])
    
    #print(houtput_list)
    
    import load
    if fnml is None:
        fnml = glob(base + '*.nml')
        if len(fnml) == 1:
            fnml = fnml[0]
            fnml = base + fnml
        else:
            fnml=input("Couldn't find any .nml file in the directory, \n "
                        "tell me which one to load: \n")
        
    nn = load.info.Nml(fnml)
    #nn.load(base + fnml)
    
    print(rsbase)
    nouts =[]
    for fn in houtput_list:
        import tree
        halo, aexp = tree.rshalo.read_halo_all(fn, sort=True, return_aexp=True, rs_dir=rsbase)
        new_aexp, nout_sim = aexp2nout(nn, aexp)
        nouts.append(nout_sim)
        print(nout_sim)
    
        if dump_pickle:
            with open(dir_out + '/halos_' + str(nout_sim).zfill(3) + '.pickle', 'wb') as fhalo:
                pickle.dump(halo, fhalo, protocol = pickle.HIGHEST_PROTOCOL)
    
        # Change .ascii file names.
        filelist = glob(fn.split('0.ascii')[0] + '*.*')
        for fn in filelist:
            former, latter = fn.split(sep='halos_')
            latter_parts = latter.split('.')
            latter_parts[0] = str(nout_sim)
            fn_new = new_dir + 'halos_' + '.'.join(latter_parts)
            copyfile(fn, fn_new)
    
    # rename out_xxx.list        
    o_lists = glob(rsbase + 'out_*.list')
    for olist in o_lists:
        with open(olist, 'r') as f:
            f.readline()        
            line = f.readline()
            aexp = float(line.split()[-1])
        new_aexp, nout_sim = aexp2nout(nn, aexp)
        former, latter = olist.split(sep='out_')
        latter_parts = latter.split('.')
        latter_parts[0] = str(nout_sim)
        fn_new = new_dir + 'out_' + '.'.join(latter_parts)
        copyfile(olist, fn_new)
    
        
#if __name__ == "__main__":
#    run()