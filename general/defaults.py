# -*- coding: utf-8 -*-
"""
Created on Mon May  9 01:23:07 2016

@author: hoseung
"""

class Default():
    
    def __init__(self):
        self.dir_snapshot = 'snapshots/'
        self.dir_halo = 'halo/'
        self.dir_galaxy = 'GalaxyMaker/'
        self.dir_halo_tree = self.dir_halo + 'Trees/'
        self.dir_galaxy_tree = self.dir_galaxy + 'Trees/'
    
    def tree_path(self, is_gal=True):
        if is_gal:
            return self.dir_galaxy + 'Trees/'
        else:
            return self.dir_halo + 'Trees/'        
        
    def ext_tree_pickle(self, is_gal=True):
        if is_gal:
            return self.dir_galaxy_tree + "extended_tree.pickle"
        else:
            return self.dir_halo_tree + "extended_tree.pickle"
        
        
