def gm2code(arr, info):
    return (arr / info.pboxsize + 0.5)

from load.rd_GM import Gal

# Tested in good/29176/
nout = 187
idgal = 90

gal = Gal(nout, idgal, load=False)
info = gal.info
# load=True by default, and 
# star = "gm", cell="gm", dm="gm".
# to set cell, dm = "raw", keep Gal() from auto-loading files.
gal.load(cell="raw", dm="raw") 

# Now gal.star, which is from GalaxyMaker dump is in gm unit
# center : (0,0,0) == (0.5,0.5,0.5) in the code unit.
# range : -0.5 * info.pboxsize ~ +0.5 * info.pboxsize
#          ~= -142 ~ +142 (at z = 0)
# and cell and dm are in the code unit.

# Unit conversion.
# move center to the galaxy center.
gal.star['x'] -= gal.header['xg'][0]
gal.star['y'] -= gal.header['xg'][1]
gal.star['z'] -= gal.header['xg'][2]

xc,yc,zc = gm2code(gal.header['xg'], info)
gal.dm['x'] -= xc
gal.dm['y'] -= yc
gal.dm['z'] -= zc

gal.cell['x'] -= xc
gal.cell['y'] -= yc
gal.cell['z'] -= zc

# rescale
gal.star['x'] *= 1e3
gal.star['y'] *= 1e3
gal.star['z'] *= 1e3

gal.dm['x'] *= info.boxtokpc
gal.dm['y'] *= info.boxtokpc
gal.dm['z'] *= info.boxtokpc

gal.cell['x'] *= info.boxtokpc
gal.cell['y'] *= info.boxtokpc
gal.cell['z'] *= info.boxtokpc

import matplotlib.pyplot as plt

fig, ax = plt.subplots(3)
ax[0].hist2d(gal.star['x'], gal.star['y'], bins=[100,100])
ax[1].hist2d(gal.dm['x'], gal.dm['y'], bins=[100,100])
ax[2].hist2d(gal.cell['x'], gal.cell['y'], bins=[100,100])
plt.show()

