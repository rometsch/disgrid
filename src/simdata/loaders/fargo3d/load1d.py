""" Functions to load 1d data from fargo3d output files.
"""
import numpy as np
import astropy.units as u

from simdata import grid
from simdata.loaders import interface


class FieldLoader1d(interface.FieldLoader):
    def load_time(self, n):
        if n is None:
            rv = self.loader.fine_output_times
        else:
            rv = self.loader.get_fine_output_time(n)
        return rv

    def load_data(self, n):
        unit = self.info["unit"]
        # + (1 if "interfaces" in self.info and "r" in self.info["interfaces"] else 0)
        Nr = self.loader.Nr
        # + (1 if "interfaces" in self.info and "phi" in self.info["interfaces"] else 0)
        Nphi = self.loader.Nphi
        if self.info["directions"] == ["r"]:
            N = Nr
        elif self.info["directions"] == ["phi"]:
            N = Nphi
        else:
            raise ValueError(
                "Trying to construct 1d field but direction is not given. Info = '{}'"
                .format(self.info))
        datafile = self.info["datafile"]
        rv = np.fromfile(datafile, count=N, offset=n * N * 8) * unit
        return rv

    def load_grid(self, n):
        r_i = np.genfromtxt(self.loader.data_dir + "/domain_y.dat"
                            )[3:-3] * self.loader.units["length"]
        # account for Fargo3d not writing out last radial interface
        if "interfaces" in self.info and "r" in self.info["interfaces"]:
            r_i = r_i[:-1]
        phi_i = np.genfromtxt(self.loader.data_dir +
                              "/domain_x.dat") * u.Unit("rad")
        active_interfaces = self.info[
            "interfaces"] if "interfaces" in self.info else []
        kwargs = {}
        for d in ["r", "phi"]:
            if d in self.info["directions"]:
                kwargs[d + "_i"] = locals()[d + "_i"]
        kwargs["active_interfaces"] = active_interfaces
        g = grid.PolarGrid(**kwargs)
        return g