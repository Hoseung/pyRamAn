# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 15:09:57 2015

@author: hoseung
"""

def static_array2d(shape, ):
    """
    array = static_array2d((100,100)))
    
    option... for fortran array
    """
    import multiprocessing
    import ctypes
    import numpy as np
    shared_array_base = multiprocessing.Array(ctypes.c_double, shape[0]*shape[1])
    # Some handles to the memory? 
    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
    # reshape => add convenience index
    return shared_array.reshape(shape[0], shape[1])