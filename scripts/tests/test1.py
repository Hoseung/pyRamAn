# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 19:45:48 2015

@author: hoseung
"""

import matplotlib.pyplot as plt
import pickle
import numpy as np

## fucntions 
def load_pickle(fname):
    with open(fname, 'rb') as f:
        return pickle.load(f)
        


wdir = '/home/hoseung/Work/data/05427'
        
cat = load_pickle(wdir + '/catalog/' + 'catalog129.pickle')

print(np.sort(cat['final_gal']))



#%%

a = [  6.08075503e+00,   2.36284615e+01,   4.21824637e+01,   8.78737116e+00,
       1.75935024e+01,   1.19991010e+02,   1.33418155e+01,   8.84852940e+00,
       7.10284997e+00,   9.17885670e+00,   3.80163984e+00,   1.28776925e+00,
       2.69120522e-01,   3.39767997e-04,   1.45614856e-04,   4.85382852e-05,
       0.00000000e+00,   2.42691426e-05,   2.42691426e-05,   0.00000000e+00,
       0.00000000e+00,   0.00000000e+00,   2.42691426e-05,   0.00000000e+00,
       4.85382852e-05]
   
   
plt.plot(a)
plt.show()