# -*- coding: utf-8 -*-

"""


From tree.fits file, select halos by mass cut.


Created on Wed Jan 21 11:45:42 2015

@author: hoseung
"""

from astropy.io import fits
from astropy.table import Table
work_dir = "/home/hoseung/Work/data/"
data = fits.getdata(work_dir + "DMO/tree/tree_eccen_v2.fits", 1)
t = Table(data)