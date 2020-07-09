""" Functions to load time data from fargo3d output files.
"""
import os
import re
import numpy as np
import astropy.units as u

from . import loadparams


def loadCoarseOutputTimes(dataDir, unit):
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