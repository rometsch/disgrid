""" Functions to load 2d data output files.
"""
import numpy as np
import astropy.units as u

from ... import grid
from .. import interface

class FieldLoader2d(interface.FieldLoader):
    def load_time(self, n, *args, **kwargs):
        if n is None:
            rv = self.loader.output_times
        else:
            if "stride" in kwargs.keys():
                n /= kwargs["stride"]
                n = int(n)
            rv = self.loader.get_time(n)
        return rv

    def load_data(self, n):
        unit = self.info["unit"]
        file_path = self.loader.data_dir + "/" + self.info["file"].format(n)
        rv = np.genfromtxt(file_path) * unit
        return rv

    def load_grid(self, n):
        xi = np.linspace(-self.loader.L/2, self.loader.L/2, self.loader.Nx+1)
        yi = np.linspace(-self.loader.L/2, self.loader.L/2, self.loader.Nx+1)

        xi = xi*self.loader.units["length"]
        yi = yi*self.loader.units["length"]
        try:
            active_interfaces = self.info["interfaces"]
        except KeyError:
            active_interfaces = []

        g = grid.CartesianGrid(x_i=xi,
                           y_i=yi,
                           active_interfaces=active_interfaces)
        return g
