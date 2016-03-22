
# coding: utf-8
def write_info(amr):
    import fortranformat as ff
    #nout = amr.nout
    
    aexp = amr.aexp
    h0 = amr.h0
    rhoc = 1.88e-29
    boxlen = 1.0
    
    f = open("info_" + str(nout).zfill(5) + ".txt", 'w')

    for name, val in zip(["ncpu", "ndim", "levelmin", "levelmax", "ngridmax", "nstep_coarse"],
                         [amr.ncpu, amr.ndim, levelmin, amr.nlevelmax, amr.ngridmax, amr.nstep_coarse]):
        f.write("{:<12s}={:11d} \n".format(name, val))
    f.write("\n")
    
    lineformat = ff.FortranRecordWriter('(1E23.15)')
    
    scale_d = amr.Om * rhoc * h0**2 / aexp**3
    scale_t = aexp**2 / (h0*1e5/3.08e24)
    scale_l = aexp* amr.boxlen * 3.08e24/(h0)
    
    for name, val in zip(["boxlen", "time", "aexp", "H0", "omega_m", "omega_l", "omega_k", "omega_b",
                         "unit_l", "unit_d", "unit_t"],
                         [boxlen, amr.t, aexp, h0, amr.Om, amr.Ol, amr.Ok, amr.Ob, scale_l, scale_d, scale_t]):
        
        f.write("{:<12s}=".format(name) + lineformat.write([val])+"\n")
    f.write("\n")
    
    f.write("ordering type=" + ah.ordering[0].decode("UTF-8"))
    
    
    f.write("\n   DOMAIN   ind_min                 ind_max  \n")
    for i in range(amr.ncpu):
        f.write("{:8d} ".format(i+1) + lineformat.write([amr.bound_key[i]]) + " "
                                   + lineformat.write([amr.bound_key[i+1]]) + "\n")
    
    f.close()

# In[1]:

wdir = "./"
import load
#info = load.info.Info()
#info.setup(nout=187, base=wdir)
nouts =[354]
for nout in nouts:
    s = load.sim.Sim(nout=nout, base=wdir)#, setup=True)
    ah = load.amr.AmrHeader()
    snout = str(nout).zfill(5)
    famr = open(wdir + "snapshots/output_"+snout+"/amr_"+snout+".out00001", 'rb')
    ah._read_amr_header(famr)
    levelmin = 7 # From other info file

    write_info(ah)

        


