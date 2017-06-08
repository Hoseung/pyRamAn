import numpy as np


def get_vmax(gal):
    if gal.vmap is None:
        get_vmap(gal)
