import numpy as np
import matplotlib
#matplotlib.use("Qt5Agg")
from matplotlib import rcParams
import matplotlib.pyplot as plt
import logging # This is a Python standard package.
import pts.simulation as sm
import pts.utils as ut
import pts.band as bnd
import pts.visual as vis

from glob import glob
from PIL import Image

import pyram

from multiprocessing import Pool

def run(gid_nout):
    name = "SDSS_RGB_LQ"
    gid, nout = gid_nout
    prefix = f"g{gid}_{nout}"
    repo = f"./{gid:05d}/{nout:05d}/"#faceon_redshift_"
    # load values from JP's value

    sim = sm.createSimulation(outDirPath=repo, prefix=prefix)
    skifile = sim.parameters()

    print("Generating mock images... gid/nout=", gid,nout)

    #for inst in sim.instruments():
    inst = sim.instruments()[0]
    totalfluxpath = inst. outFilePaths(fileType="total.fits")[0]
    datacube = sm.loadFits(totalfluxpath) # loadFits & getFitsAxes return 'astropy' quantities with units attached.
    x,y,wavelengths = sm.getFitsAxes(totalfluxpath)

    if bands is not None:
        # 0.6, 0.75, 1.2
        # 0.5, 0.8, 1.5 at z=0.3
        contributions = [ (bands[0], 0.5, 0, 0), (bands[1], 0, 0.8, 0), (bands[2], 0, 0, 1.5) ]
        # Could loop over sims
        # Make RGB images of ALL instruments.
        _, _ = vis.makeConvolvedRGBImages(sim,
                                        contributions=contributions,
                                        fileType=output_type,
                                        name=name,
                                        stretch='log',
                                        decades=decs[0],
                                        fmin=fmin_log1,
                                        fmax=fmax_log1)
        _, _ = vis.makeConvolvedRGBImages(sim,
                                        contributions=contributions,
                                        fileType=output_type,
                                        name=name,
                                        stretch='log',
                                        decades=decs[1],
                                        fmin=fmin_log2,
                                        fmax=fmax_log2)

simdir = '/scratchb01/hoseung/'
wdir = './'
rgb_done = glob(wdir+"faceon_35/g*_*_faceon_total_SDSS_RGB_LQ_log_dec3.5.png")
has_fits = glob(wdir+"?????/?????/*face_on_total.fits")

output_type = "total"
quality = ["intermediate", "high"][0]
colors = "SDSS_I,SDSS_G,SDSS_U" #MASS_2MASS_H,2MASS_2MASS_J,2MASS_2MASS_KS"
nthreads = 1

do_all=True
if do_all:
    fits_fns =[(int(fn.split('/')[-3]), int(fn.split('/')[-2])) for fn in has_fits]
    fits_fns.sort()

    rgb_done_fns=[]
    for fn in glob(wdir+"faceon_35/g*_*_faceon_total_SDSS_RGB_LQ_log_dec3.5.png"):
        ss = fn.split('faceon_35/g')[1].split("_")
        rgb_done_fns.append((int(ss[0]),int(ss[1])))
    rgb_done_fns.sort()
    
    fits_only = [(gid,nout) for (gid,nout) in fits_fns if (gid,nout) not in rgb_done_fns]
else:
    fits_only = [(27, 874), (1, 874)]
print(fits_only)

wavelengths = None
if wavelengths is not None:
    tuples = { name: wavelengths << sm.unit("micron") }
    for sim in sims:
        vis.makeRGBImages(sim, wavelengthTuples=tuples, fileType=type)


decs=[3.5,4.0]
# Keep flux range the same throughout snapshots.z`
fmin_ash, fmax_ash =  2e-3, 2e2
fmin_log1, fmax_log1 = 1e-2, 1e2
fmin_log2, fmax_log2 = 1e-3, 1e2

segments = colors.split(',')
if len(segments) != 3:
    raise ut.UserError("colors argument must have three comma-separated segments")
try: bands = [ bnd.BroadBand(segment) for segment in segments ]
except ValueError: bands = None

with Pool(nthreads) as p:
    p.map(run, fits_only)

