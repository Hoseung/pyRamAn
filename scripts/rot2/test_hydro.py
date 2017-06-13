import load

s = load.sim.Sim(nout=782, base='/data52/Horizon-AGN/OUTPUT_DIR/')
s.set_ranges(ranges=[[0.4,0.41]]*3)
s.add_hydro()
