# -*- coding: utf-8 -*-
"""
Mass histogram of selected clusters
together with halo mass function.

Created on Wed Jun 24 17:30:09 2015

@author: hoseung
"""

import csv
import numpy as np

targets = '/home/hoseung/Copy/Results/zoominTarget/targets.txt'

with open(targets, 'r') as f:
    a = f.readlines()
    ntargets = len(a) -2

id = np.zeros(ntargets)
mvir=np.zeros(ntargets)
rvir=np.zeros(ntargets)
x=np.zeros(ntargets)
y=np.zeros(ntargets)
z=np.zeros(ntargets)
data=[]
with open(targets, 'r') as f:
    f.readline()
    for i, line in enumerate(f.readlines()):
        ll = line.split()
        if len(ll) > 1:
            id[i] = ll[0]
            mvir[i] = ll[1]
            rvir[i] = ll[2]
            x[i] = ll[3]
            y[i] = ll[4]
            z[i] = ll[5]

#        data.append((line.split())

#%%
import matplotlib.pyplot as plt
nbin = 4
n, bins, patches = plt.hist(np.log10(mvir), nbin, 
                            range=(13.5,15.5), 
                            align='mid',
                            normed=False, 
                            facecolor='green', 
                            alpha=0.75)
plt.locator_params(nbins=5)
plt.xlabel('Cluster mass [$M_{\odot}$ ]', fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.ylabel('# Clusters', fontsize = 20)
plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)
xx=[13.5, 14.0, 14.5, 15.0, 15.5]
yy=[30, 25, 20, 15, 10]
plt.plot(xx,yy)
plt.show()