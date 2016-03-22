# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 15:38:32 2015

@author: hoseung
"""

import numpy as np
import csv 

class Allrefine():
    def __init__(self, hid, data, aexp):
        self.id = hid
        self.x = data['x_refine']
        self.y = data['y_refine']
        self.z = data['z_refine']
        self.r = data['r_refine']
        self.aexp = aexp

class TypedWriter:
    """ 
    A CSV writer which will write rows to CSV file "f",
    which uses "fieldformats" to format fields.
    """
    def __init__(self, f, fieldnames, fieldformats, **kwds):
        self.writer = csv.DictWriter(f, fieldnames, **kwds) 
        self.formats = fieldformats

    def writerow(self, row):
        self.writer.writerow(dict((k, self.formats[k] % v) for k, v in row.iteritems()))

    def writerows(self, rows):
        for row in rows:
            self.writerow(row)

def plot_dist3d(x, y, z, rr, save=None, show=False):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
#    ax = fig.add_subplot(1,1,1, projection = '3d')
    
# line plot
    N = len(x)
    for i in range(N-1):
        ax.plot(x[i:i+2], y[i:i+2], z[i:i+2], color=plt.cm.jet(255*i/N))
    ax.scatter(x,y,z, rr)
# cx, cy, cz makes the data 2D-like.
# ... HOW???
    cx = np.ones_like(x) * ax.get_xlim3d()[0]
    cy = np.ones_like(x) * ax.get_ylim3d()[1]
    cz = np.ones_like(z) * ax.get_zlim3d()[0]

    for i in range(N-1):
        ax.plot(x[i:i+2], y[i:i+2], cz[i:i+2], color=plt.cm.jet(255*i/N))
        ax.plot(x[i:i+2], cy[i:i+2], z[i:i+2], color=plt.cm.jet(255*i/N))
        ax.plot(cx[i:i+2], y[i:i+2], z[i:i+2], color=plt.cm.jet(255*i/N))

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.locator_params(nbins = 4)

    if not show:
        plt.ioff()

    if save is not None:
        fn_out = save + '3d.png'
        plt.savefig(fn_out, dpi=300)

    if show:
        plt.show()

    plt.close()


#%%
#
# interpolate into NOUT (read all aouts) and interpolate linearly.
# aexp_new = np.arange(0.,1.0,0.005) # (start,stop,step) 
# or 
#rscale=raw_input('Multiply radius?')
#if not isinstance(rscale,int) :
#    rscale = 1.0
#print "R_scale =",rscale
work_dir = './'
rscale = 1.0
nnouts = 201
nzoomlevel = 2
music=False
n_more_level = 8
aexp_new = np.linspace(0.0,1.0,nnouts)
#aexp_org = np.arange(0,1,0.0125) + 0.0125  # nout = 80
nout_ini = 1
nout_fi = 11
import pickle 
f_rf = work_dir + 'refine_params.pickle'
with open(f_rf, "rb") as f:
    rf = pickle.load(f)

print(rf.x, rf.y, rf.z)

    
#%%
param_dir = work_dir + 'refine_params/'

ok = rf.x[0] > 0.00001

# recent version of ramses read refine parameters of each levels.

x_refine_new = np.zeros((len(aexp_new),nzoomlevel + n_more_level), dtype='float')
y_refine_new = np.zeros((len(aexp_new),nzoomlevel + n_more_level), dtype='float')
z_refine_new = np.zeros((len(aexp_new),nzoomlevel + n_more_level), dtype='float')
r_refine_new = np.zeros((len(aexp_new),nzoomlevel + n_more_level), dtype='float')
m_refine_new = np.zeros((len(aexp_new),nzoomlevel + n_more_level), dtype='float')
coollev = np.array([1]*len(aexp_new), dtype='int')

for hind in range(len(rf.id)):

    if music:
        x_off = rf.x[hind][-1] - 0.5 # last position goes to 0.5
        y_off = rf.y[hind][-1] - 0.5
        z_off = rf.z[hind][-1] - 0.5
        rf.x[hind] = rf.x[hind] - x_off
        rf.y[hind] = rf.y[hind] - y_off
        rf.z[hind] = rf.z[hind] - z_off
    
    x_refine_new[:,nzoomlevel] = np.interp(aexp_new,rf.aexp,rf.x[hind][ok])
    y_refine_new[:,nzoomlevel] = np.interp(aexp_new,rf.aexp,rf.y[hind][ok])
    z_refine_new[:,nzoomlevel] = np.interp(aexp_new,rf.aexp,rf.z[hind][ok])
    r_refine_new[:,nzoomlevel] = np.interp(aexp_new,rf.aexp,rf.r[hind][ok])*rscale
  
    galname = str(rf.id[hind]).zfill(6)
    
    for i in range(nzoomlevel):
        for j in range(len(aexp_new)):
            x_refine_new[j,i] = x_refine_new[j,nzoomlevel]
            y_refine_new[j,i] = y_refine_new[j,nzoomlevel]
            z_refine_new[j,i] = z_refine_new[j,nzoomlevel]
            r_refine_new[j,i] = r_refine_new[j,nzoomlevel] + (nzoomlevel - i) * 0.01
            m_refine_new[j,i] = 0.0
    
    for i in range(nzoomlevel, nzoomlevel + n_more_level):
        for j in range(len(aexp_new)):
            x_refine_new[j,i] = x_refine_new[j,nzoomlevel]
            y_refine_new[j,i] = y_refine_new[j,nzoomlevel]
            z_refine_new[j,i] = z_refine_new[j,nzoomlevel]
            r_refine_new[j,i] = r_refine_new[j,nzoomlevel]
            m_refine_new[j,i] = 8.0
#    plot_dist3d(x_refine_new, y_refine_new, z_refine_new, r_refine_new,
#                save=param_dir + galname, show=False)
#    plot_dist3d(rf.x[hind][ok],rf.y[hind][ok],rf.z[hind][ok],rf.r[hind][ok],
#                save=param_dir + galname + "aa")    
    x_refine_new = x_refine_new[::-1]
    y_refine_new = y_refine_new[::-1]
    z_refine_new = z_refine_new[::-1]
    r_refine_new = r_refine_new[::-1]
    m_refine_new = m_refine_new[::-1]
    with open(param_dir + galname + 'refine_params.txt', 'w') as f:
        writer = csv.writer(f)
        print("nnout", end="\n", file=f)
        print(nnouts, end="\n", file=f)     
        print("nzoomlevel", end="\n", file=f)
        print(nzoomlevel, end="\n", file=f)
        print("AEXP", end="\n", file=f)
        for i in aexp_new[:-1]:
            print("{:4f}".format(i), end=", ", file=f)
        print("%.4f" % aexp_new[-1], end="\n", file=f)
        print("X_REFINE", end="\n", file=f)
        for item in x_refine_new:
            for ii in item:
                print("{:5f} ".format(ii), end=' ', file=f)
            print(" ", file=f)

        print("Y_REFINE", end="\n", file=f)
        for item in y_refine_new:
            for ii in item:
                print("{:5f} ".format(ii), end=' ', file=f)
            print(" ", file=f)        

        print("Z_REFINE", end="\n", file=f)
        for item in z_refine_new:
            for ii in item:
                print("{:5f} ".format(ii), end=' ', file=f)
            print(" ", file=f)

        print("R_REFINE", end="\n", file=f)
        for item in r_refine_new:
            for ii in item:
                print("{:5f} ".format(ii), end=' ', file=f)
            print(" ", file=f)

        print("M_REFINE", end="\n", file=f)
        for item in m_refine_new:
            for ii in item:
                print("{:5f} ".format(ii), end=' ', file=f)
            print(" ", file=f)
            
        print("COOLLEV", end="\n", file=f)
        for ii in coollev:
            print(ii, file=f)
