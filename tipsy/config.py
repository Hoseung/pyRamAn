import numpy as np
class Changa_config():
    def __init__(self):
        self.gas_properties = ["iord", "Metalsdot", "OxMassFrac","OxMassFracdot", 
                 "coolontime", "HeI", "HeII", "HI", "rung"]
        self.star_properties = ["iord", "rung", "igasorder", "massform", 
                  "Metalsdot", "OxMassFrac","OxMassFracdot"]
        self.dm_properties = ["iord", "rung"]
        self.dtype_aux = np.dtype([("iord", "<i4"), 
                                  ("rung", "<i4"),
                                  ("igasorder", "<i4"), 
                                  ("massform", "<f4"),
                                  ("Megalsdot", "<f4"),
                                  ("OxMassFrac", "<f4"),
                                  ("OxMassFracdot", "<f4"),
                                  ("coolontime", "<f4"),
                                  ("HeI", "<f4"),
                                  ("HeII", "<f4"),
                                  ("HI", "<f4")])


        self.dtype_gas = np.dtype({  'm': (('<f4', 1), 0), 
                         'pos': (('<f4', (3,)), 4), 
                           'x': (('<f4', 1), 4), 
                           'y': (('<f4', 1), 8), 
                           'z': (('<f4', 1), 12),
                         'vel': (('<f4', (3,)), 16),
                          'vx': (('<f4', 1), 16),
                          'vy': (('<f4', 1), 20),
                          'vz': (('<f4', 1), 24),
                         'rho': (('<f4', 1), 28),
                        'temp': (('<f4', 1), 32),
                         'eps': (('<f4', 1), 36),
                       'metal': (('<f4', 1), 40),
                         'phi': (('<f4', 1), 44)})

        self.dtype_dm = np.dtype({  'm': (('<f4', 1), 0), 
                         'pos': (('<f4', (3,)), 4), 
                           'x': (('<f4', 1), 4), 
                           'y': (('<f4', 1), 8), 
                           'z': (('<f4', 1), 12),
                         'vel': (('<f4', (3,)), 16),
                          'vx': (('<f4', 1), 16),
                          'vy': (('<f4', 1), 20),
                          'vz': (('<f4', 1), 24),
                         'eps': (('<f4', 1), 28),
                         'phi': (('<f4', 1), 32)})

        self.dtype_star = np.dtype({  'm': (('<f4', 1), 0), 
                         'pos': (('<f4', (3,)), 4), 
                           'x': (('<f4', 1), 4), 
                           'y': (('<f4', 1), 8), 
                           'z': (('<f4', 1), 12),
                         'vel': (('<f4', (3,)), 16),
                          'vx': (('<f4', 1), 16),
                          'vy': (('<f4', 1), 20),
                          'vz': (('<f4', 1), 24),
                       'metal': (('<f4', 1), 28),
                       'tform': (('<f4', 1), 32),
                         'eps': (('<f4', 1), 36),
                         'phi': (('<f4', 1), 40)})
