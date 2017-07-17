import numpy as np
from collections import Iterable
import utils.match as mtc

"""
    NOTE
    ----
    Instead of using if ... elif clause, you can use decorator to implement multiple dispatch.
    But That's more appropriate for fundamental operations.
    See an example below. (Ffrom http://matthewrocklin.com/blog/work/2014/02/25/Multiple-Dispatch)

    >>> from multipledispatch import dispatch

    >>> @dispatch(int, int)
    ... def add(x, y):
    ...     return x + y

    >>> @dispatch(object, object)
    ... def add(x, y):
    ...     return "%s + %s" % (x, y)

"""


class Nnza():
    def __init__(self, fname="./nout_nstep_zred_aexp.txt"):
        self.nnza = np.genfromtxt(fname,
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])
    def step2out(self, nstep):
        if isinstance(nstep, int):
            return int(self.nnza["nout"][np.where(self.nnza["nstep"] == nstep)[0]])
        elif not isinstance(nstep, str) and isinstance(nstep, Iterable):
            return self.nnza["nout"][mtc.match_list_ind(self.nnza["nstep"], nstep)]
            

    def out2step(self, nout):
        if isinstance(nout, int):
            return int(self.nnza["nstep"][np.where(self.nnza["nout"] == nout)[0]])
        elif not isinstance(nout, str) and isinstance(nout, Iterable):
            return self.nnza["nstep"][mtc.match_list_ind(self.nnza["nout"], nout)]

