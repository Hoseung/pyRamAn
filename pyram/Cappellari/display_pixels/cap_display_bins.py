"""
NAME:
    display_bins()
    
AUTHOR:
    Michele Cappellari, University of Oxford
    cappellari_at_astro.ox.ac.uk

PURPOSE:
    This simple routine illustrates how to display a Voronoi binned map.
    
INPUTS:
    (x, y): (length npix) Coordinates of the original spaxels before binning;
    binNum: (length npix) Bin number corresponding to each (x, y) pair,
            as provided in output by the voronoi_2d_binning() routine;
    velBin: (length nbins) Quantity associated to each bin, resulting
            e.g. from the kinematic extraction from the binned spectra.
          
MODIFICATION HISTORY:
    V1.0.0: Michele Cappellari, Oxford, 15 January 2015          
    
"""

from cap_display_pixels import display_pixels

def display_bins(x, y, binNum, velBin):

    npix = len(binNum)
    if (npix != len(x)) or (npix != len(y)):
        raise ValueError('The vectors (x, y, binNum) must have the same size')
        
    f = display_pixels(x, y, velBin[binNum])
    
    return f
