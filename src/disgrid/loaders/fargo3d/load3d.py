""" Functions to load 3d data from fargo3d output files.
"""
import numpy as np
import astropy.units as u

from disgrid import grid
from disgrid.loaders import interface

from .loadbinary import load_data

class FieldLoader3d(interface.FieldLoader):
    def load_time(self, n, *args, **kwargs):
        if n is None:
            rv = self.loader.output_times
        else:
            if "stride" in kwargs.keys():
                n /= kwargs["stride"]
                n = int(n)
            rv = self.loader.get_output_time(n)
        return rv

    def load_data(self, n):
        unit = self.info["unit"]
        # + (1 if "interfaces" in self.info and "r" in self.info["interfaces"] else 0)
        Nr = self.loader.Nr
        # + (1 if "interfaces" in self.info and "phi" in self.info["interfaces"] else 0)
        Nphi = self.loader.Nphi
        Ntheta = self.loader.Ntheta
        file_path = self.loader.data_dir + "/" + self.info["pattern"].format(n)
        rv = load_data(file_path, self.info["varname"], n)
        rv = rv.reshape(Ntheta, Nr, Nphi)
        rv = np.swapaxes(rv, 0, 2)
        rv = np.swapaxes(rv, 0, 1)
        return rv *unit

    def load_grid(self, n):
        r_i = np.genfromtxt(self.loader.data_dir + "/domain_y.dat"
                            )[3:-3] * self.loader.units["length"]
        # account for Fargo3d not writing out last radial interface
        if "interfaces" in self.info and "r" in self.info["interfaces"]:
            r_i = r_i[:-1]
        phi_i = np.genfromtxt(self.loader.data_dir +
                              "/domain_x.dat") * u.Unit("rad")
        theta_i = np.genfromtxt(self.loader.data_dir + "/domain_z.dat"
                            )[3:-3] * u.Unit("rad")
        active_interfaces = self.info[
            "interfaces"] if "interfaces" in self.info else []
        g = grid.SphericalGrid(r_i=r_i,
                           phi_i=phi_i,
                           theta_i=theta_i,
                           active_interfaces=active_interfaces)
        return g
