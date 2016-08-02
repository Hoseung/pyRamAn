
# coding: utf-8

# In[5]:

import matplotlib.pyplot as plt 


class Merger():
    def __init__(self):
        pass

import pickle
mpgs = pickle.load(open(wdir + "main_prgs.pickle", "rb"))


# In[11]:

wdir = './29172/'
for gg in mpgs:
    from matplotlib.backends.backend_pdf import PdfPages
    fig, ax = plt.subplots(1, sharex=True)
    galidx = gg.data["idx"][i]
    with PdfPages(wdir + 'stellar_density_map_'+ str(galidx) + '.pdf') as pdf:

        for i, nout in enumerate(gg.nouts):
            galid = gg.data["id"][i]
            gal = Gal(nout, galid, wdir=wdir, load=False)
            gal.load(star="raw", cell="raw", dm="raw")
            gal.star["x"] = gal.star["x"] - gal.header["xg"][0]
            gal.star["y"] = gal.star["y"] - gal.header["xg"][1]
            gal.star["z"] = gal.star["z"] - gal.header["xg"][2]
            ax.hist2d(gal.star["x"], gal.star["y"], bins=[100,100])
            ax.set_title("id: {} at nout: {}".format(galid, nout))
            pdf.savefig()
            ax.clear()

        #info = gal.info
        
        #print(gal.data["xc"][0] / gal.data[""])

