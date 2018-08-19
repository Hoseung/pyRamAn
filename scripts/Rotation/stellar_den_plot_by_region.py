
# coding: utf-8

# In[5]:

import matplotlib.pyplot as plt 
import pickle
from matplotlib.backends.backend_pdf import PdfPages
from load.rd_GM import Gal

class Merger():
    def __init__(self):
        pass

import os 
wdir = './'
outdir = wdir + "stellar_density_evol/"
if not os.path.isdir(outdir):
    os.mkdir(outdir)

rscale = 1.0
mpgs = pickle.load(open(wdir + "main_prgs.pickle", "rb"))
radius = 0.001

fig, ax = plt.subplots(1, sharex=True)
for gg in mpgs:
    try:
        with PdfPages(outdir + 'stellar_density_map_'+ \
                      str(gg.ids[0]) + "_" + \
                      str(gg.idxs[0].astype(int)) + '.pdf') as pdf:

            for i, nout in enumerate(gg.nouts):
                galid = gg.ids[i]
                if gg.data["rgal"][i] > 0:
                    gal = Gal(nout, galid, base=wdir, load=False)
    #                print("nout and Galaxy size", nout, gal.header["r"])
                    gal.load(star="raw", cell="none", dm="none",
                             rscale=rscale, radius=radius)
                    gal.star["x"] = gal.star["x"] - gal.header["xg"][0]
                    gal.star["y"] = gal.star["y"] - gal.header["xg"][1]
                    gal.star["z"] = gal.star["z"] - gal.header["xg"][2]
                    ax.hist2d(gal.star["x"], gal.star["y"], bins=[100,100])
                    ax.set_title("id: {} at nout: {}".format(galid, nout))
                    ax.set_xlim([-radius,radius])
                    ax.set_ylim([-radius,radius])
                    ax.set_aspect('equal', 'datalim')
                    pdf.savefig()
                    ax.clear()
                else:
                    print("skip ---------", nout)
    except:
        pass

