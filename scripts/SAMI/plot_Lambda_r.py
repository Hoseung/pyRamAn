# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 18:19:56 2015

@author: hoseung
"""

def plot_lambda(catalog, i_early, i_late, i_bad, fn_out='./'):
    import matplotlib.pyplot as plt
    plt.ioff()
    f = plt.figure(figsize=(8,6))

    ax = f.add_subplot(111)
    #for i, val in enumerate(lambdar_arr):
    for i in i_early:
        a = np.asarray(catalog['lambda_arr'][i][:8])
        ax.plot(a, 'r-', alpha=0.5) # Red = Early
#    for i in i_late:
#        ax.plot(catalog['lambda_arr'][i], 'b-', alpha=0.3) # Red = Early
    
    #plt.xlabel() # in the unit of Reff
    ax.set_title("56 galaxies", fontsize=22) 
    ax.set_ylabel(r"$\lambda _{R}$", fontsize=22) 
    ax.set_xlabel(r'$R/R_{eff}$', fontsize=16)
    ax.set_xlim(right=7)
    ax.set_xticks([0, 3.5, 7])
    ax.set_xticklabels(["0", "0.5", "1"])
#    plt.savefig(fn_out)
#    plt.close()    
    
    
import pickle
import numpy as np

wdir = '/home/hoseung/Work/data/05427/'

nout_fi = 187

def load_pickle(fname):
    with open(fname, 'rb') as f:
        return pickle.load(f)

# Final galaxies
# Why set..?

#for inout, nout in enumerate(reversed(np.arange(nout_ini, nout_fi + 1))):
cat = load_pickle(wdir + 'catalog/' + 'catalog' + str(nout_fi) + '.pickle')

i_early = np.where(cat['mstar'] > 0)[0]
i_early = list(np.arange(len(cat)))
i_late = []
i_bad = [1, 12, 30, 50]
#bad_gals = [1496, 85, 1636, 1340]

for i in i_bad:
    i_early.remove(i)


#    if not os.path.isdir(wdir + snout + '/'):
#        os.mkdir(wdir + snout + '/')
import matplotlib.pylab as plt
plot_lambda(cat, i_early, i_late, i_bad, fn_out = wdir + "187_lambdar_disk.png")
plt.show()
#%%