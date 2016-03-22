
# coding: utf-8

# In[13]:

import os
repo = input("repository: \n")

full_path = os.path.realpath(repo)
nout_ini = int(input("first nout: \n"))
nout_fi = int(input("last nout: \n"))
import sys

if os.path.exists(full_path + 'GalaxyMaker/inputfiles_HaloMaker.dat'):
    print("Warnning, the file already exists. \n Do you want to overwrite?")
    OK = input("y/n")
    if OK != "y":
        sys.exit()
        
with open("inputfiles_HaloMaker.dat", 'w') as f:
    for nout in range(nout_ini, nout_fi + 1):
        f.write("'" + full_path + '/snapshots/output_' + str(nout).zfill(5))
        f.write("'  Ra3  1   " + str(nout) + "\n")
        
        


# In[ ]:



