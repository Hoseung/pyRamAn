# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 17:53:03 2015

Check refine param 


@author: hoseung
"""

def read_params(fname):
    f = open(fname, 'r')
    f.readline() # skip the first line
    nout=int(f.readline())
    f.readline()
    nzoomlevel=f.readline()
    f.readline()
    #aexp = float(f.readline())
    s=f.readline() # s is a buch of 'comma-separated' values
    aexp = [float(x) for x in s.split(',')]

    x_refine=[]
    y_refine=[]
    z_refine=[]
    r_refine=[]

    #x_refine
    f.readline()
    for i in range(0,nout):
        x_refine.append(float(f.readline()))
    #y_refine
    f.readline()
    for i in range(0,nout):
        y_refine.append(float(f.readline()))
    #z_refine
    f.readline()
    for i in range(0,nout):
        z_refine.append(float(f.readline()))
    f.readline()
    for i in range(0,nout):
        r_refine.append(float(f.readline()))
    f.close()
    
    return x_refine, y_refine, z_refine, r_refine
    
# First refine param
import glob
import matplotlib.pyplot as plt
for file in glob.glob("*refine*.txt"):
    idgal = file[0:5]
    print(file)
    x1, y1, z1, r1 = read_params(file)
# Second one
#x1, y1, z1, r1 = read_params("refines_01191.txt")
# parameter 1
    plt.ioff()
    fig = plt.figure()
    ax1 = fig.add_subplot(2,2,1)
    ax1.plot(x1, c='blue')
    
    ax2 = fig.add_subplot(2,2,2)
    ax2.plot(y1)
    
    ax3 = fig.add_subplot(2,2,3)
    ax3.plot(z1)
    
    ax4 = fig.add_subplot(2,2,4)
    ax4.plot(r1)
    plt.savefig(idgal+".png")

"""
##################################
#parameter 2
    ax1.plot(x2, c='red')
    ax2.plot(y2, c='red')
    ax3.plot(z2, c='red')
    ax4.plot(r2, c='red')
"""
    