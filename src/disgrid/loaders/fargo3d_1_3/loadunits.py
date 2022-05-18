""" Functions to load units from fargo3d output files.
"""
import os
import re

import numpy as np
import astropy.units as u
import astropy.constants as const

from . import loadparams


def loadUnits(dataDir):
    units = {
        "time" : 1*u.s,
        "mass" : 1*u.g,
        "length" : 1*u.cm
    }
    return units
