
# coding: utf-8

# Table 1 of the paper.
# List of clusters, 
# mass, radius, number of galaxies, (Can I measure relaxedness?)
# and so on if needed later.

# In[38]:

import load
import pickle
import tree.halomodule as hmo
import numpy as np


# load halo bricks of each cluster, 
# galaxyMaker brick (for total number of galaxies)
# galaxy catalog (for number of galaxies above the mass cut)

# In[51]:

def get_ngals(halo, gal):  
    return np.where(np.square(gal.data.x - halo.x) + 
                    np.square(gal.data.y - halo.y) +
                    np.square(gal.data.z - halo.z) < np.square(halo.rvir))[0]
    
    
import glob
base = './'

#clusters = glob.glob(base + "?????")
clusters = ['01605', '07206', '35663', '24954', '49096',\
            '05427', '05420', '29172', '29176', '10002',\
            '36415', '06098', '39990', '36413', '17891', '04466']

print(clusters)

main_halos=[]
ngals=[]
len_cat=[]
len_gal_all=[]

info = load.info.Info(base=clusters[0], nout=187, load=True)

for cluster in clusters:
    wdir = cluster + '/'
    try:
        nout = 187
        h = hmo.Halo(nout=nout, base=wdir, halofinder='HM', is_gal=False, load=True)
    except:
        nout = 170
        h = hmo.Halo(nout=nout, base=wdir, halofinder='HM', is_gal=False, load=True)
    hg = hmo.Halo(nout=nout, base=wdir, halofinder='HM', is_gal=True, load=True)
    cat = pickle.load(open(wdir + 'easy_final/catalog' + str(nout) + '.pickle', 'rb'))
    clu = h.data[h.data.np.argmax()]
    main_halos.append(clu)
    gals_in_cluster = get_ngals(clu,hg)
    ngals.append(len(gals_in_cluster))
    len_cat.append(len(cat))
    len_gal_all.append(len(hg.data))
    
    print(cluster, clu.mvir/1e13, clu.rvir * info.pboxsize, clu.nsub, len(cat))


# In[60]:


f = open(base + 'cluster_list.txt', 'w')
f.write("ID      Mvir     Rvir     # member    #large gals  # all gals\n")
for a,b,c,d,e in zip(clusters, main_halos, ngals, len_cat, len_gal_all,):
    f.write("{}  {}  {}  {}  {}  {}\n".format(a, b.mvir/1e13, b.rvir*info.pboxsize, c, d,e))
    print(a, b.mvir/1e13, b.rvir*info.pboxsize, c, d, e)


# In[ ]:



