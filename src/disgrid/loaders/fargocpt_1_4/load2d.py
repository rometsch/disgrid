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
            rv = self.loader.snapshot_time(n)
        return rv

    def load_data(self, n):
        unit = u.Unit(self.info["unit"])
        Nr = self.info["Nrad"]
        Naz = self.info["Nazi"]
        filepath = f"snapshots/{n}/" + self.info["filename"]
        filepath = self.loader.datadir_path(filepath)
        rv = np.fromfile(filepath).reshape(Nr, Naz) * unit
        return rv

    def load_grid(self, n):
        active_interfaces = []
        if self.info["on_radial_interface"]:
            active_interfaces.append("r")
        if self.name == "vazi":
            active_interfaces.append("phi")
        g = grid.PolarGrid(r_i=self.loader.r_i,
                           phi_i=self.loader.phi_i,
                           active_interfaces=active_interfaces)
        return g
