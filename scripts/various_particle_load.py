# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 15:42:34 2015

To test various particle load options

@author: hoseung
"""

# 1. DMO
# 1.1 mass / id / velocity / ref


# 2. General
# 2.1 dm:  mass / id / velocity / ref
# 2.2 dm & star  mass / id / velocity / ref
# 2.3 dm & sink  mass / id / velocity / ref
# 2.4 star & sink  mass / id / velocity / ref

# Think about speed! 


base = '/home/hoseung/Work/data/C01605_mid/'
nout = 124
ptype = "dm"
frefine= 'refine_params.txt'
fnml = 'cosmo_200.nml'
scale = 0.01

import load

snout = str(nout).zfill(3)

s = load.sim.Sim(nout, base)
ptype=["dm mass id", "star  star vel mAss"]
s.add_part(ptypes = ptype)
s.set_ranges([[0.5,0.55],[0.5,0.55],[0.5,0.55]])

s.part.load()