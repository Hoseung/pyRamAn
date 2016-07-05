import load

s = load.sim.Sim(base='./', nout=187)
s.add_part(ptypes=["star id pos vel", "dm id pos"])
s.add_hydro()


import tree.halomodule as hmo

hals = hmo.Halo(nout=187, is_gal=False)
gals = hmo.Halo(nout=187, is_gal=True)
