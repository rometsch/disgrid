""" Functions to load 1d data from fargo3d output files.
"""
import os

import astropy.units as u
import numpy as np
from ... import grid
from .. import interface


def load1dRadial(dataFile, unit, lengthunit=None):
    data = np.fromfile(dataFile, dtype=float)
    if 'torque_planet' in dataFile:
        v = data[1::2] * unit
        r = data[::2]
    else:
        v = data[1::4] * unit
        r = data[::4]
    if lengthunit is None:
        return v
    else:
        r = r * lengthunit
        return (r, v)


def load1dMassFlow(n, dataDir):
    # parse datafile info
    Nr = 0
    unit = 1
    with open(os.path.join(dataDir, "gasMassFlow1D.info"), "r") as infofile:
        for line in infofile:
            line = line.strip()
            parts = [s.strip() for s in line.split("=")]
            if parts[0] == "Nr":
                Nr = int(parts[1])
            elif parts[0] == "unit":
                unit = u.Unit(parts[1])
    # get data
    Nsamples = 100
    data = np.fromfile(os.path.join(dataDir, "gasMassFlow1D.dat"),
                       count=Nsamples * Nr,
                       offset=Nr * n * 8)
    data = data.reshape([int(len(data) / Nr), Nr])
    data = np.average(data, axis=0)
    if n == 0:
        data = np.zeros(len(data))
    rv = data * unit
    return rv


def parse_1d_info_file(fname):
    stem = os.path.basename(os.path.abspath(fname))[:-7]
    pattern = "{}1D{}.dat".format(stem, "{}")
    with open(fname, "r") as infofile:
        for line in infofile:
            line = line.strip()
            parts = [s.strip() for s in line.split("=")]
            if parts[0] == "Nr":
                Nr = int(parts[1])
            elif parts[0] == "unit":
                unit = parts[1]
    return {"unit": unit, "Nr": Nr, "pattern": pattern}


class FieldLoader1d(interface.FieldLoader):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fileinfo = parse_1d_info_file(
            self.loader.datadir_path(self.info["infofile"])
        )

    def load_time(self, n):
        if n is None:
            rv = self.loader.output_times
        else:
            rv = self.loader.snapshot_time(n)
        return rv

    def load_data(self, n):
        unit = u.Unit(self.fileinfo["unit"])
        datafile = self.loader.datadir_path(self.fileinfo["pattern"].format(n))
        rv = load1dRadial(datafile, unit)
        return rv

    def load_grid(self, n):
        if self.fileinfo["Nr"] == len(self.loader.r_i):
            active_interfaces = ["r"]
        else:
            active_interfaces = []
        g = grid.PolarGrid(r_i=self.loader.r_i,
                           active_interfaces=active_interfaces)
        return g


class FieldLoader1dTorq(interface.FieldLoader):
    def load_time(self, n):
        if n is None:
            rv = self.loader.output_times
        else:
            rv = self.loader.snapshot_time(n)
        return rv

    def load_data(self, n):
        datafile = self.loader.datadir_path(self.info["pattern"].format(n))
        data = np.fromfile(datafile, dtype=float)
        rv = data[1::2] * u.Unit("cm^2 g / s^2")
        return rv

    def load_grid(self, n):
        g = grid.PolarGrid(r_i=self.loader.r_i)
        return g


class FieldLoader1dMassFlow(interface.FieldLoader):
    def load_time(self, n):
        if n is None:
            rv = self.loader.fine_output_times
        else:
            rv = self.loader.get_fine_output_time(n)
        return rv

    def load_data(self, n):
        rv = load1dMassFlow(n, self.loader.data_dir)
        return rv

    def load_grid(self, n):
        r_i = np.genfromtxt(self.loader.data_dir +
                            "/used_rad.dat") * u.Units(self.loader.units["length"])
        g = grid.PolarGrid(r_i=r_i, active_interfaces=["r"])
        return g
