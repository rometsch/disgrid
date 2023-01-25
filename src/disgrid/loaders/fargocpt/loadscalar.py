""" Functions to scalar time data from fargo3d output files.
"""
import os
from functools import lru_cache

import astropy.units as u
import numpy as np
from ... import scalar


@lru_cache(100)
def load_text_data_variables(filepath):
    # load all variable definitions from a text file
    # which contains the variable names and colums in its header.
    # each variable is indicated by
    # "#variable: {column number} | {variable name} | {unit}"
    found_variables = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line[0] != "#":
                break
            identifier = "#variable:"
            if line[:len(identifier)] == identifier:
                col, name, unitstr = [
                    s.strip() for s in line[len(identifier):].split("|")
                ]
                found_variables[name] = (col, unitstr)
    return found_variables

def remove_duplicates(vecs):
    """Remove duplicate entries from time series.
    
    The first array in the vecs list is assumed to be the time/index array.
    All duplicate entries in this array are identified and removed in the time/index array
    and all corresponding entries in the remaining arrays in vecs are also removed.
    
    This proocess is repeated recursively to remove lines which might be repeated multiple times.

    Parameters
    ----------
    vecs: list of arrays
        A list of the arrays to remove duplicates in.

    Returns
    -------
    vecs: list of arrays
    """
    x = vecs[0]
    identical = x[1:] == x[:-1]
    if np.sum(identical) == 0:
        return vecs
    inds = np.append([True], np.logical_not(identical))
    new_vecs = [v[inds] for v in vecs]
    return remove_duplicates(new_vecs)

@lru_cache(100)
def load_data(filepath):
    return np.genfromtxt(filepath).T


def load_text_data_file(filepath, varname, Nmax=np.inf):
    # get data
    variables = load_text_data_variables(filepath)
    col = int(variables[varname][0])
    unit_str = variables[varname][1]
    unit_str = unit_str.replace("1/s", "s-1")
    unit = u.Unit(unit_str)
    file_data = load_data(filepath)
    data = file_data[col] * unit
    time_col = int(variables["physical time"][0])
    time = file_data[time_col]
    N = min(len(data), len(time))
    data = data[:N]
    time = time[:N]
    time, data = remove_duplicates([time, data])
    if data.isscalar:
        data = u.quantity.Quantity([data])
    N = min(len(data), Nmax)
    data = data[:N]
    return data


class ScalarLoader:
    def __init__(self, name, datafile, loader, *args, **kwargs):
        self.loader = loader
        self.datafile = datafile
        self.name = name

    def __call__(self):
        time = self.load_time()
        data = self.load_data()
        f = scalar.Scalar(time, data, name=self.name)
        return f

    def load_data(self):
        datafile = self.loader.datadir_path(self.datafile)
        rv = load_text_data_file(datafile,
                                 self.name,
                                 Nmax=len(self.loader.fine_output_times))
        return rv

    def load_time(self):
        datafile = self.loader.datadir_path(self.datafile)
        rv = load_text_data_file(datafile,
                                 "physical time",
                                 Nmax=len(self.loader.fine_output_times))
        return rv
