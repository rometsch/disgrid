""" Functions to scalar time data from fargo3d output files.
"""
import os
from functools import lru_cache
import logging

import astropy.units as u
import numpy as np
from ... import scalar

@lru_cache(20)
def _load_text_data_variables(filepath, timestamp):
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
    logging.debug(found_variables)
    return found_variables

def load_text_data_variables(filepath):
    return _load_text_data_variables(filepath, os.path.getmtime(filepath))


@lru_cache(20)
def _load_data(filepath, timestamp):
    return np.genfromtxt(filepath).T

def load_data(filepath):
    return _load_data(filepath, os.path.getmtime(filepath))


def load_text_data_file(filepath, varname, Nmax=np.inf):
    # get data
    variables = load_text_data_variables(filepath)
    col = int(variables[varname][0])
    unit_str = variables[varname][1]
    unit_str = unit_str.replace("1/s", "s-1")
    unit = u.Unit(unit_str)
    file_data = load_data(filepath)
    data = file_data[col] * unit
    time_col = int(variables["time"][0])
    time = file_data[time_col]
    try:
        N = min(len(data), len(time))
        N = min(len(data), Nmax)
        data = data[:N]
        time = time[:N]
    except TypeError:
        N = 1
    if data.isscalar:
        data = u.quantity.Quantity([data])
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
        rv = self.loader.fine_output_times
        return rv
