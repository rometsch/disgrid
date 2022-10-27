""" Functions to load time data from output files.
"""
import os
import numpy as np


def load_times(dataDir, unit):
    """ Load the time array.

    Parameters
    ----------
    dataDir: str
        Path to the data directory.
    unit: astropy.units.Unit, timelike
        Code unit of the time.
        
    Returns
    -------
    astropy.units.Quantity array
        Time corresponding to the full outputs.
    """
    times_file = os.path.join(dataDir, "times.dat")
    times = np.genfromtxt(times_file) * unit
    return times