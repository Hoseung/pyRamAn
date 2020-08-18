import numpy as np
class Config():
    """
    All string variables except for path are stored in lowercsae.
    """
    def __init__(self, 
                cosmo=True,
                sim_type=None,
                dmo=False,
                longint=False,
                sink_path='./',
                snap_dir='./',
                dtype_c=None
                ):
        self.cosmo = cosmo
        self.sim_type = sim_type
        self.dmo = dmo
        self.longint = longint
        self.sink_path=sink_path
        self.snap_dir=snap_dir
        self.dtype_c=dtype_c

    @property
    def cosmo(self):
        return self._cosmo
    @cosmo.setter
    def cosmo(self, cosmo):
        self._cosmo=cosmo

    @property
    def sim_type(self):
        return self._sim_type
    @sim_type.setter
    def sim_type(self, sim_type):
        try:
            self._sim_type=sim_type.lower()
        except:
            pass

    @property
    def dmo(self):
        return self._dmo
    @dmo.setter
    def dmo(self, dmo):
        self._dmo=dmo

    @property
    def longint(self):
        return self._longint
    @longint.setter
    def longint(self, longint):
        self._longint=longint

    @property
    def sink_path(self):
        return self._sink_path
    @sink_path.setter
    def sink_path(self, sink_path):
        self._sink_path=sink_path

    @property
    def dtype_c(self):
        return self._dtype_c
    @dtype_c.setter
    def dtype_c(self, dtype_c):
        self._dtype_c=np.dtype(dtype_c)

    @property
    def snap_dir(self):
        """
        This is equivalent with Simbase._fnbase
        Need to be merged
        """
        return self._snap_dir
    @snap_dir.setter
    def snap_dir(self, snap_dir):
        self._snap_dir=snap_dir

    

def empty_config():
    conf = {'cosmo':None,
            'sim_type':None,
            'dmo':None,
            'longint':None,
            'sink_path':'./',
            }
    return conf
