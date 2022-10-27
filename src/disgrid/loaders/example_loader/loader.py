import copy
import os
import glob
import re

import astropy.constants as const
import astropy.units as u
import numpy as np
from disgrid import field, fluid, grid, particles
from disgrid.loaders import interface

from . import defs
from . import load2d
from . import loadparams
from . import loadtime
from . import loadunits
from . import loadscalar

code_info = ("example", "0.0", "description")

def identify(path):
    """ Identifies a directory to be a example simulation dir.

    Parameters
    ----------
    path: str
        Path of directory to be checked.
    
    Returns
    -------
    bool
        True if dir is a example code dir, False if not.
    """
    idfile_path = os.path.join(path, "codeid.txt")
    if os.path.exists(idfile_path):
        with open(idfile_path, "r") as infile:
            content = infile.read().strip()
            if content == "example code":
                return True
    else:
        return False


class Loader(interface.Interface):

    code_info = code_info

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data_dir = self.path
        self.output_times = np.array([])

    def scout(self):
        self.get_parameters()
        self.get_units()
        self.apply_units()
        self.load_times()
        self.get_domain_size()
        self.get_planets()
        self.get_fluids()
        self.get_scalars()
        self.get_fields()
        # self.get_nbodysystems()
        # self.register_alias()

    def get_parameters(self):
        self.parameters = loadparams.load_params(self.path)

    def get_nbodysystems(self):
        pass

    def apply_units(self):
        for vardict in [defs.planet_vars_scalar, defs.vars_2d, defs.vars_scalar]:
            for info in vardict.values():
                if "unit" in info:
                    info["unit"] = u.Unit(info["unit"])
                else:
                    info["unit"] = loadunits.get_unit_from_powers(info["unitpowers"],
                                                    self.units)

    def get_planets(self):
        # create planets
        self.planets.append(particles.Planet("planet", 0))
        planet = self.planets[0]
        # add variables to planets:
        for varname in defs.planet_vars_scalar:
            info = defs.planet_vars_scalar[varname]
            datafile = os.path.join(self.data_dir, info["file"])
            loader = loadscalar.ScalarLoader(varname, datafile, info, self)
            planet.register_variable(varname, loader)

    def get_fluids(self):
        fluid_names = os.listdir(os.path.join(self.data_dir,"fluids"))
        for name in fluid_names:
            self.fluids[name] = fluid.Fluid(name)

    def get_fields(self):
        self.get_fields_2d()

    def get_fields_2d(self):
        for fluidname, fl in self.fluids.items():
            for varname, info in defs.vars_2d.items():
                info = copy.deepcopy(info)
                info["file"] = f"fluids/{fluidname}/" + info["file"]
                info["varname"] = varname
                fieldLoader = load2d.FieldLoader2d(varname, info, self)
                fl.register_variable(varname, "2d", fieldLoader)

    def get_scalars(self):
        for fluid_name in self.fluids:
            fl = self.fluids[fluid_name]
            fluid_data_dir = os.path.join(self.data_dir, "fluids", fluid_name)
            for varname, info in defs.vars_scalar.items():
                datafile = os.path.join(fluid_data_dir, info["file"])
                if os.path.exists(datafile):
                    fl.register_variable(
                        varname, "scalar",
                        loadscalar.ScalarLoader(varname, datafile, info, self))

    def get_domain_size(self):
        self.Nx = self.parameters["Nx"]
        self.Ny = self.parameters["Ny"]
        self.L = self.parameters["L"]

    def load_times(self):
        self.output_times = loadtime.load_times(self.data_dir, self.units["time"])

    def get_time(self, n):
        return self.output_times[n]

    def get_units(self):
        self.units = loadunits.load_units(self.data_dir)


