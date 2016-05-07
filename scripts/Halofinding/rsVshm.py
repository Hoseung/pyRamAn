# -*- coding: utf-8 -*-
"""
Test resemblance of HM and RS trees

Created on Sun Jun 14 00:33:52 2015

@author: hoseung
"""

# Load HM
from astropy.io import fits
from astropy.table import Table
import tree
wdir = '/home/hoseung/Work/data/AGN2/'
data = fits.getdata(wdir + "halo/TMtree.fits", 1)
hmt = Table(data)

#%%

idgals = [5232, 5495, 5543, 6061, 5491, 6191]
for idgal in idgals:
    prg_treex, prg_tree = tree.tmtree.get_main_prg(hmt, idgal, nout_ini=122, nout_fi=0)    
    print(prg_treex)

#%%
# Load RS


# 3D plot
tree.treeplots(hmt, thisgal, save=save_dir)




#%%
all_final_halo = tru.final_halo_list(hmt)
## ID list

#%%
# mass-produce plots of halo properties.
quantities=["sam_mvir","mvir","rvir","rs","vrms","vmax"
    ,"jx","jy","jz","spin","m200b","m200c","m500c","m2500c"
    ,"xoff","voff","btoc","ctoa","ax","ay","az"]

normalizer=np.array([1e-11,1e-11,1,1,1,1,1,1,1e-11
        ,1,1e-11,1e-11,1e-11,1e-11,1,1,1,1,1,1,1])

for i, hal in enumerate(all_final_halo[0:10]):
    print(i, hal)
    tree = tru.get_main_prg_tree(hmt, hal)

    fn_save = str(hal) + 'halo_all.pdf'
#    trp.plot_all(tree, hal, save=True, out_dir=work_dir, fn_save=fn_save,
#                 nrows=4, ncols=math.ceil(len(quantities)/4),
#                 quantities=quantities, normalizer=normalizer)
#    trp.plot_all_multiPDF(tree, hal, out_dir=work_dir + 'RS_trees/', fn_save=fn_save,
#                 nrows=2, ncols=2,
#                 quantities=quantities, normalizer=normalizer)
    trp.trajectory3D(tree, hal, save=work_dir + 'RS_trees/')