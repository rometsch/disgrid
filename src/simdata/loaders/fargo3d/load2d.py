""" Functions to load 2d data from fargo3d output files.
"""
import os
import numpy as np
import astropy.units as u

from simdata import grid
from simdata.loaders import interface


class FieldLoader2d(interface.FieldLoader):
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
        # + (1 if "interfaces" in self.info and "phi" in self.info["interfaces"] else 0)

        if self.loader.Ntheta > 1:
            shape = (self.loader.Ntheta, self.loader.Nr)
        else:
            shape = (self.loader.Nr, self.loader.Nphi)

        data = np.fromfile(os.path.join(self.loader.data_dir,
                            self.info["pattern"].format(n))).reshape(*shape) * unit

        return data

    def load_grid(self, n):
        r_i = np.genfromtxt(self.loader.data_dir + "/domain_y.dat"
                            )[3:-3] * self.loader.units["length"]
        # account for Fargo3d not writing out last radial interface
        if "interfaces" in self.info and "r" in self.info["interfaces"]:
            r_i = r_i[:-1]

        active_interfaces = self.info[
            "interfaces"] if "interfaces" in self.info else []

        if self.loader.Ntheta > 1:
            theta_i = np.genfromtxt(os.path.join(self.loader.data_dir, "domain_z.dat"))[3:-3] * u.Unit("rad")
            g = grid.SphericalGrid(r_i=r_i,
                                    theta_i=theta_i,
                                    active_interfaces=active_interfaces)
        else:
            phi_i = np.genfromtxt(self.loader.data_dir +
                                  "/domain_x.dat") * u.Unit("rad")
            g = grid.PolarGrid(r_i=r_i,
                            phi_i=phi_i,
                            active_interfaces=active_interfaces)

        return g
