# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 01:12:12 2015

@author: hoseung
"""

def static_array2d(shape):
    import multiprocessing
    import ctypes
    import numpy as np
    shared_array_base = multiprocessing.Array(ctypes.c_double, shape[0]*shape[1], lock=False)
#    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
#    shared_array_base = mp.RawArray(ctypes.c_double, shape[0]*shape[1])

    # Some handles to the memory? 
    shared_array = np.ctypeslib.as_array(shared_array_base)
    # reshape => add convenience index
    return shared_array.reshape(shape[0], shape[1])
    




def mk_gal(halodata, out_q, info, i, dm, star, cell,
            save=False, rscale=0.3, verbose=True):
#            , def_param=(dm, star, cell)):
    from galaxymodule import galaxy

#    print(halodata['id'])
    print("{}-th **********************************************".format(i))


    
    
    




import time
t1 = time.time()

m = mp.Manager()
out_q = m.Queue()

nh = len(h.data)

pool = mp.Pool(processes=ncore)
for i in range(nh):
    print(i,"-th started")
    pool.apply(mk_gal, args=(h.data[i], out_q, s.info, i, dm, star, cell))

time.sleep(2)
print("----------Done---------")
print(" Took", time.time() - t1, "Seconds")