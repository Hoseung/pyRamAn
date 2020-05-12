
# coding: utf-8
def write_info(amr):
    #import fortranformat as ff
    #nout = amr.nout
    
    aexp = amr.aexp
    h0 = amr.h0 * 1e-2
    rhoc = 1.88e-29
    boxlen = 1.0
    
    f = open("info_" + str(nout).zfill(5) + ".txt", 'w')

    for name, val in zip(["ncpu", "ndim", "levelmin", "levelmax", "ngridmax", "nstep_coarse"],
                         [amr.ncpu, amr.ndim, levelmin, amr.nlevelmax, amr.ngridmax, amr.nstep_coarse]):
        f.write("{:<12s}={:11d} \n".format(name, val))
    f.write("\n")
    
    #lineformat = ff.FortranRecordWriter('(1E23.15)')

    scale_d = amr.Om * rhoc * h0**2 / aexp**3 
    scale_t = aexp**2 / (h0*1e5/3.08e24)
    scale_l = aexp* amr.boxlen * 3.08e24/(h0) 
    
    for name, val in zip(["boxlen", "time", "aexp", "H0", "omega_m", "omega_l", "omega_k", "omega_b",
                         "unit_l", "unit_d", "unit_t"],
                         [boxlen, amr.t, aexp, h0, amr.Om, amr.Ol, amr.Ok, amr.Ob, scale_l, scale_d, scale_t]):
        
        f.write("{:<12s}={:.15E} \n".format(name,val))
    f.write("\n")
    
    f.write("ordering type=" + ah.ordering[0].decode("UTF-8"))
    
    
    f.write("\n   DOMAIN   ind_min                 ind_max  \n")
    for i in range(amr.ncpu):
        f.write("{:8d}   {:.15E}   {:.15E}\n".format(i+1, amr.bound_key[i],amr.bound_key[i+1]))
    
    f.close()

"""
This can generate 'header' of info. 
But it is not trivial to read 128-bit floating point (QUADHILBERT) numbers from binary bits in Python. 
Instead, I used a fortran program to read amr.00001 and output hilbert keys in the info format.
"""

wdir = "./"
from pyram import load
nouts = range(113, 120)
for nout in nouts:
    ah = load.sim.AmrHeader()
    snout = str(nout).zfill(5)
    ah._read_amr_header(open(wdir + "output_"+snout+"/amr_"+snout+".out00001", 'rb'), skip_header=False)
    levelmin = 8 # From other info file

    write_info(ah)

        


