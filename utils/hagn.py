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
        nnza = np.genfromtxt(fname,
                     dtype=[("nout", int),
                            ("nstep", int),
                            ("zred", float),
                            ("aexp", float)])

        self.nnza = np.zeros(len(nnza),dtype=[("nout", int),
                                    ("nstep", int),
                                    ("zred", float),
                                    ("aexp", float),
                                    ("lbt", float)])
        self.nnza["nout"] = nnza["nout"]
        self.nnza["nstep"] = nnza["nstep"]
        self.nnza["zred"] = nnza["zred"]
        self.nnza["aexp"] = nnza["aexp"]

    def a2b(self, vals, field_a, field_b):
        if not (field_a in self.nnza.dtype.names  and  field_b in self.nnza.dtype.names):
            raise ValueError("One or both of {} {} not found".format(field_a, field_b))
        if isinstance(vals, np.integer):
            return int(self.nnza[field_b][np.where(self.nnza[field_a] == vals)[0]])
        elif not isinstance(vals, str) and isinstance(vals, Iterable):
            return self.nnza[field_b][mtc.match_list_ind(self.nnza[field_a], vals)]
        else:
            NotImplementedError()

    def step2out(self, nstep):
        """
        Deprecated, Use a2b instead.
        """
        if isinstance(nstep, np.integer):
            return int(self.nnza["nout"][np.where(self.nnza["nstep"] == nstep)[0]])
        elif not isinstance(nstep, str) and isinstance(nstep, Iterable):
            return self.nnza["nout"][mtc.match_list_ind(self.nnza["nstep"], nstep)]
        else:
            NotImplementedError()


    def out2step(self, nout):
        """
        Deprecated, Use a2b instead.
        """
        if isinstance(nout, np.integer):
            return int(self.nnza["nstep"][np.where(self.nnza["nout"] == nout)[0]])
        elif not isinstance(nout, str) and isinstance(nout, Iterable):
            return self.nnza["nstep"][mtc.match_list_ind(self.nnza["nout"], nout)]
        else:
            NotImplementedError()
