import numpy as np

class Nnza():
    def __init__(self, fname="./nout_nstep_zred_aexp.txt"):
        self.nnza = np.genfromtxt(fname,
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])
    def step2out(self, nstep):
        return int(self.nnza["nout"][np.where(self.nnza["nstep"] == nstep)[0]])

    def out2step(self, nout):
        return int(self.nnza["nstep"][np.where(self.nnza["nout"] == nout)[0]])


