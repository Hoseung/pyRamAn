# -*- coding: utf-8 -*-
"""
Created on Fri May 29 02:43:32 2015


Pickle Rockstar ASCII outputs.

@author: hoseung
"""
import tree
from os.path import isdir
from os import path
from os import mkdir
from os import rename
import pickle
from glob import glob

base = '/home/hoseung/Work/data/AGN2/'

# directory that pickled halos go
dir_out = path.normpath(base + 'rhalo/halos_py')
if not isdir(dir_out):
    mkdir(dir_out)

rsbase = base + 'rhalo/rockstar_halos/'
print("go")

#%%
# directory that re-named files go
mkdir(rsbase + 'new')
fout = open(rsbase + 'new/datasets.txt', 'w')
with open(rsbase + 'datasets.txt', 'r') as fdataset:
    temp = fdataset.readline()
    fout.write(temp)
    for s in fdataset.readlines():
        sa, sb = s.split()
        nout_sim = int(sa.split('/info_')[1].split('.txt')[0])
        nout_rs = int(sb)

        # write the modified line to the new file.
        fout.write(sa + '\t' + str(nout_sim)+ '\n')

        # load halo
        fn = rsbase + 'halos_' + str(nout_rs) + '.0'
        halo = tree.rshalo.read_halo_all(fn, sort=True)
        # and pickle halo.
        with open(dir_out + '/halos_' + str(nout_sim).zfill(3) + '.pickle' , 'wb') as fhalo:
            pickle.dump(halo, fhalo, protocol = pickle.HIGHEST_PROTOCOL)
            
        # Change .ascii file names.
        filelist = glob(rsbase + 'halos_' + str(nout_rs) + '.*.*')
        for fn in filelist:
            former, latter = fn.split(sep='halos_')
            latter_parts = latter.split('.')
            latter_parts[0] = str(nout_sim)
            fn_new = former + 'new/halos_' + '.'.join(latter_parts)
            rename(fn, fn_new)

fout.close()