""" Functions to scalar time data from fargo3d output files.
"""
import numpy as np
from disgrid import scalar 

class ScalarLoader:
    def __init__(self, name, datafile, info, loader, *args, **kwargs):
        self.loader = loader
        self.datafile = datafile
        self.info = info
        self.name = name
        self.units = loader.units

    def __call__(self):
        time = self.load_time()
        data = self.load_data()
        f = scalar.Scalar(time, data, name=self.name)
        return f

    def load_data(self):
        col = self.info["datacol"]
        unit = self.info["unit"]
        rv = np.genfromtxt(self.datafile, usecols=int(col)) * unit
        return rv

    def load_time(self):
        col = self.info["timecol"]
        unit = self.units["time"]
        rv = np.genfromtxt(self.datafile, usecols=int(col)) * unit
        return rv



