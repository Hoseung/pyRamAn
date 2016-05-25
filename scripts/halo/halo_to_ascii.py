
# coding: utf-8

# In[7]:

def halo_to_ascii(hh, fname, quick=False):
    with open(fname, 'w') as f:
        if quick:
            for line in hh.data.dtype.descr:
                f.write(line[0] + " ")
            f.write("\n")
            hh.data.tofile(f, sep=" \n")
        else:
            for line in hh.data.dtype.descr:
                f.write(line[0] + " ")
            f.write("\n")

            for line in hh.data:
                for ll in line:
                    f.write(str(ll) + ' ')
                f.write('\n')


# In[8]:

import load
import tree.halomodule as hmo

wdir = "./"
for nout in range(20, 188):
    hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', is_gal=False)
    fname = wdir + 'halo_brick_' + str(nout).zfill(3) + '.txt'
    halo_to_ascii(hh, fname, quick=False)

