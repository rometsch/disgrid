""" Functions to load 2d data from fargocpt output files.
"""
import numpy as np
from astropy import units as u
from ... import grid
from .. import interface


class FieldLoader2d(interface.FieldLoader):
    def load_time(self, n):
        if n is None:
            rv = self.loader.output_times
        else:
            rv = self.loader.get_output_time(n)
        return rv

    def load_data(self, n):
        if "unit" in self.info:
            unit = u.Unit(self.info["unit"])
        else:
            unit = 1
            for baseunit, power in self.info["unitpowers"].items():
                unit = unit * u.Unit(baseunit)**power
        Nr = self.loader.Nr + (1 if "interfaces" in self.info
                               and "r" in self.info["interfaces"] else 0)
        Nphi = self.loader.Nphi  #+ (1 if "interfaces" in self.info and "phi" in self.info["interfaces"] else 0)
        filepath = self.loader.filepath(self.info["pattern"].format(n))
        rv = np.fromfile(filepath).reshape(Nr, Nphi) * unit
        return rv

    def load_grid(self, n):
        active_interfaces = self.info[
            "interfaces"] if "interfaces" in self.info else []
        g = grid.PolarGrid(r_i=self.loader.r_i,
                           phi_i=self.loader.phi_i,
                           active_interfaces=active_interfaces)
        return g
