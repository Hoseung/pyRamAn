class Merger():
    pass

import os
import shutil
import pickle
wdir = "./"
cdir = "easy_final/"
mpgs = pickle.load(open(wdir + "main_prgs.pickle", "rb"))

for gallist in mpgs:
    galdir = wdir + cdir + str(gallist.ids[0]) + "_" + str(gallist.idxs[0])
    if not os.path.isdir(galdir):
        os.mkdir(galdir)
    for nout, idgal in zip(gallist.nouts, gallist.ids):
        fn = str(nout) + "_" + str(idgal) + ".png"
        #flist.append(fn)
        try:
            os.link(wdir + cdir + "galaxy_plot/" + fn, galdir + "/" + fn)
        except:
            pass
            # Phantom galaxies don't have figures. 
    
