""" Functions to load units from fargo3d output files.
"""
import os
import re

import numpy as np
import astropy.units as u
import astropy.constants as const

from . import loadparams


def loadUnits(dataDir):
    # load data units
    first_summary = os.path.join(
        dataDir, loadparams.find_first_summary(dataDir))
    if os.path.exists(first_summary):
        ptrn = r"COMPILATION OPTION SECTION:\n==============================\n.*\-DCGS.*\nGhost"
        with open(first_summary, 'r') as infile:
            if re.search(ptrn, infile.read()):
                # have cgs units
                units = {
                    "mass": u.g,
                    "time": u.s,
                    "length": u.cm,
                    "temperature": u.K
                }
                return units

        # Try to extract unit normalisation from summary
        units = {}
        units["mass"] = 1.0
        units["time"] = 1.0
        units["length"] = 1.0

        with open(first_summary, 'r') as infile:
            ptrn = r"(?<=R0 = \()\d+.\d+"
            m = re.search(ptrn, infile.read())
            if m:
                units["length"] *= float(m.group()[0])

        with open(first_summary, 'r') as infile:
            ptrn = r"(?<=MSTAR = \()\d+.\d+"
            m = re.search(ptrn, infile.read())
            if m:
                units["mass"] *= float(m.group())

        with open(first_summary, 'r') as infile:
            ptrn = r"STEFANK =.*\*pow\(\(\d+\.\d+\)\/\((\d+\.*\d*\*\d+\.\d*\w+\d*)\),-0\.5"
            m = re.search(ptrn, infile.read())
            if m:
                components = m.group(1).split('*')
                units["length"] *= float(components[0]) * float(
                    components[1]) * u.cm

        with open(first_summary, 'r') as infile:
            ptrn = r"STEFANK =.*\*pow*\(\(\d+\.\d*\)\/(\d+\.\d*\w+\d*),-1\.5\)"
            m = re.search(ptrn, infile.read())
            if m:
                components = m.group(1).split('*')
                units["mass"] *= float(m.group(1)) * u.g

        units["time"] = (np.sqrt(units["length"]**3 /
                                 (const.G.cgs * units["mass"]))).to(u.s)

    return units
