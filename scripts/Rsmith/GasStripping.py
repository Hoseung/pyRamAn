# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 22:55:17 2015

Calculate the amount of gas stripped. 

@author: hoseung
"""


# 1. Prepare galaxy

## define halo 

# ->> Done. 
import halo
h = halo.aaa

## Galaxy as baryons within 0.5 Rvir

class Smbh():
    def __init__(self, sink):
        self.cloud = sink
        self.mbh = m * ncloud
    
    def load_cloud(self, cloud):
        self.set_com(x, y, z)

    def set_com(x, y, z):
        """
        Calculate the center of mass of cloud particles 
        as the mean position of cloud particles.
        
        notes
        -----
            This is not precisely the center of mass. 
            How does RAMSES define C.O.M? 
            Cloud particles are all at the same mass. No need to consider.
        
        """

class Galaxy():
    def __init__(self):
        pass
        
    def add_dm(self, dm):
        self.dm = dm
        pass

    def add_star(self, star):
        self.star = star
        pass

    def add_sink(self, sink):
        self.smbh = Smbh(sink)



# 2. measure gas stripping

#  distinguish hot and cold gas (10^4K)

# M_cold = sum(cell.rho[temperature < 1e4] * cell.dx[temperature < 1e4]**3)
# M_hot = sum(cell.rho[temperature > 1e4] * cell.dx[temperature > 1e4]**3)