# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 14:46:30 2015

@author: hoseung

HM halo util

"""

def print_halo(halo):
    for names in halo.dtype.names:
        print(names, halo[names][0])  #, data[names].shape)

def norm_halo(halo, info):
    # To do: Check if it is already normalized or not
    halo['p'][0][0] = halo['p'][0][0] / info.pboxsize + 0.5
    halo['p'][0][1] = halo['p'][0][1] / info.pboxsize + 0.5
    halo['p'][0][2] = halo['p'][0][2] / info.pboxsize + 0.5
    halo['r'][0] = halo['r'][0] / info.pboxsize
    halo['rvir'][0] = halo['rvir'][0] / info.pboxsize
    halo['m'][0] = halo['m'][0] * 1e11

def load_data(nout, work_dir=None, normalize=True):
    snout = str(nout).zfill(3)
    try:
        from scipy.io.idl import readsav
        data = readsav(work_dir + 'halo/halo' + snout + '.sav')['h']
        if normalize is True:
            import load
            info = load.info.Info()
            info.setup(nout=nout, base=work_dir)
            info.read_info()
            norm_halo(data, info)

    except:
        print("Cannot specify a file to read")
        print("trying to read {0}".format(work_dir + 'halo/halo' + snout + '.sav'))
        print("+++++++++++++++++++++++")
    # Readsav returns a recarray.


    return data
