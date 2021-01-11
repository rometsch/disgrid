""" Functions to load time data from fargo3d output files.
"""
import os
import re
import numpy as np
import astropy.units as u

from simdata.monotonize import monotonize
from . import loadparams


def loadCoarseOutputTimes(dataDir, unit):
    """ Load the time array corresponding to the full output steps.
    
    If there is a planet in the simulation, we can extract the information
    from its data file. Otherwise we need to look into all summary files.
    
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
    planet_file = os.path.join(dataDir, "planet0.dat")
    if os.path.exists(planet_file):
        return loadCoarseOutputTimesFromPlanet(planet_file, unit)
    else:
        return loadCoarseOutputTimesFromSummary(dataDir, unit)
        
def loadCoarseOutputTimesFromPlanet(planet_file, unit):
    """ Load the time array from the planet file.
    
    Parameters
    ----------
    planet_file: str
        Path to the planet output file.
    unit: astropy.units.Unit, timelike
        Code unit of the time.
        
    Returns
    -------
    astropy.units.Quantity array
        Time corresponding to the full outputs.
    """
    time = np.genfromtxt(planet_file, usecols=8)
    try:
        inds = monotonize(time)
        time = time[inds]
    except IndexError:
        time = np.array([time])
    return time * unit

def loadCoarseOutputTimesFromSummary(dataDir, unit):
    # search all summary.dat files for the time
    outputTimes = []
    pattern = re.compile(r'OUTPUT [0-9]* at simulation time ([0-9\.eE+-]*)')
    for f in sorted([f for f in os.listdir(dataDir) if 'summary' in f],
                    key=lambda x: int(x[7:-4])):
        with open(os.path.join(dataDir, f), 'r') as infile:
            datastr = infile.read()
            matches = re.findall(pattern, datastr)
            try:
                outputTimes.append(float(matches[0]))
            except ValueError:
                break
    times = np.array(outputTimes)
    return times * unit


def loadFineOutputTimes(dataDir, unit):
    numbers = loadparams.find_summary_numbers(dataDir)
    times = np.array([])
    for n in numbers:
        params = loadparams.getParamsFromNthSummary(dataDir, n)
        dt = params["dt"]
        Ninterm = params["ninterm"]
        offset = 0 if len(times) == 0 else times[-1]
        new_times = np.arange(1, Ninterm + 1) * dt + offset
        times = np.append(times, new_times)
    times = times * unit
    return times