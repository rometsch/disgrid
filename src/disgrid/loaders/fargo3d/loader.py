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
from . import load1d
from . import load2d
from . import load3d
from . import loadparams
from . import loadtime
from . import loadunits
from . import loadgrid
from . import loadscalar

code_info = ("fargo3d", "2.0", "public")


def identify(path):
    """ Identifies a directory to be a fargo3d simulation dir.

    Parameters
    ----------
    path: str
        Path of directory to be checked.
    
    Returns
    -------
    bool
        True if dir is a fargo3d dir, False if not.
    """
    try:
        get_data_dir(path)
        return True
    except FileNotFoundError:
        return False


def var_in_files(varpattern, files):
    """ Search a list of filenames for an output file name pattern.

    Paramters
    ---------
    varpattern: str
        Pattern to identify a quantity output file as defined in the definition file.

    Returns
    -------
    bool
        True if found, False if not.
    """
    p = re.compile(varpattern.replace(".", r"\.").format(r"\d+"))
    for f in files:
        if re.match(p, f):
            return True
    return False


def get_data_dir(path):
    rv = None
    ptrn = re.compile(r"summary\d+.dat")
    for glob_pattern in ["*", "*/*", "*/*/*", "*/*/*/*", "**"]:
        for f in glob.iglob(os.path.join(path, glob_pattern)):
            m = re.search(ptrn, f)
            if m:
                rv = os.path.dirname(f)
                break
        if rv is not None:
            break
    if rv is None:
        raise FileNotFoundError(
            r"Could not find identifier file 'summary\d+.dat' in any subfolder of '{}'"
            .format(path))
    return rv

def parse_power_def(unit_def, dim):
    exponent = 0
    if isinstance(unit_def, str):
        unit_def = [unit_def]
    try:
        iter(unit_def)
    except TypeError:
        unit_def = [unit_def]
    for p in unit_def:
        if p == "maxdim":
            p = dim
        elif p =="-maxdim":
            p = -dim
        exponent += p
    return exponent

def get_unit_from_powers(unitpowers, units, dim):
    unit = 1.0
    for k, p in unitpowers.items():
        expo = parse_power_def(p, dim)
        unit = unit * units[k]**expo
    return unit


class Loader(interface.Interface):

    code_info = code_info

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data_dir = get_data_dir(self.path)
        self.output_times = np.array([])
        self.fine_output_times = np.array([])

    def scout(self):
        self.get_parameters()
        self.get_domain_size()
        self.calc_dimension()
        self.get_units()
        self.apply_units()
        self.load_times()
        self.get_planets()
        self.get_fluids()
        self.get_fields()
        self.get_scalars()
        self.get_nbodysystems()
        self.register_alias()

    def get_parameters(self):
        self.parameters = loadparams.getParamsFromNthSummary(
            self.data_dir, loadparams.find_first_summary_number(self.data_dir))

    def apply_units(self):
        for vardict in [defs.planet_vars_scalar, defs.vars_maxdim, defs.vars_1d, defs.vars_scalar]:
            for info in vardict.values():
                if "unit" in info:
                    info["unit"] = u.Unit(info["unit"])
                else:
                    info["unit"] = get_unit_from_powers(info["unitpowers"],
                                                    self.units, self.dim)

    def register_alias(self):
        for particlegroup in self.particlegroups:
            particlegroup.alias.register_dict(defs.alias_particle)
        for planet in self.planets:
            planet.alias.register_dict(defs.alias_particle)
        for f in self.fluids.values():
            f.alias.register_dict(defs.alias_fields)
            f.alias.register_dict(defs.alias_reduced)

    def get_nbodysystems(self):
        pass

    def get_planets(self):
        planet_ids = []
        p = re.compile(r"bigplanet(\d).dat")
        for s in os.listdir(self.data_dir):
            m = re.match(p, s)
            if m:
                planet_ids.append(m.groups()[0])
        planet_ids.sort()
        # create planets
        self.planets = []
        for pid in planet_ids:
            self.planets.append(particles.Planet(str(pid), pid))
        # add variables to planets
        for pid, planet in zip(planet_ids, self.planets):
            for varname in defs.planet_vars_scalar:
                info = defs.planet_vars_scalar[varname]
                datafile = os.path.join(self.data_dir,
                                        info["file"].format(pid))
                loader = loadscalar.ScalarLoader(varname, datafile, info, self)
                planet.register_variable(varname, loader)

    def get_fluids(self):
        ptrn = re.compile(r"output(.*)\.dat")
        fluid_names = [
            m.groups()[0]
            for m in (re.search(ptrn, f) for f in os.listdir(self.data_dir))
            if m is not None
        ]
        for name in fluid_names:
            self.fluids[name] = fluid.Fluid(name)

    def get_fields(self):
        with open(os.path.join(self.data_dir, "outputgas.dat"), "r") as in_file:
            self.mpiio_vars = [line.strip().split()[1][3:] for line in in_file]
        if self.dim == 3:
            self.get_fields_multidim("3d", load3d.FieldLoader3d)
        elif self.dim == 2:
            self.get_fields_multidim("2d", load3d.FieldLoader2d)
        self.get_fields_1d()

    def get_fields_multidim(self, dim, loader_class):
        files = os.listdir(self.data_dir)
        for fluidname in self.fluids.keys():
            for varname, info in defs.vars_maxdim.items():
                info_formatted = copy.deepcopy(info)
                info_formatted["pattern"] = info_formatted["pattern"].format(
                    fluidname, "{}")
                info_formatted["varname"] = varname
                if var_in_files(info_formatted["pattern"], files):
                    fieldLoader = loader_class(varname, info_formatted, self)
                    self.fluids[fluidname].register_variable(
                        varname, dim, fieldLoader)
            for varname in self.mpiio_vars:
                pretty_name = defs.mpiio_vars[varname]
                info_formatted = copy.deepcopy(
                    defs.vars_maxdim[pretty_name])
                info_formatted["pattern"] = "{}_{}.mpiio".format(
                    fluidname, "{}")
                info_formatted["varname"] = varname
                fieldLoader = loader_class(pretty_name, info_formatted, self)
                self.fluids[fluidname].register_variable(
                    pretty_name, dim, fieldLoader)

    def get_fields_1d(self):
        for fluid_name in self.fluids:
            fl = self.fluids[fluid_name]
            monitor_dir = os.path.join(self.data_dir, "monitor", fluid_name)
            for name_pattern, info in defs.vars_1d.items():
                for n in range(len(self.planets)):
                    try:
                        filename = info["pattern"].format(n)
                    except IndexError:
                        filename = info["pattern"].format(fluid_name, n)

                    datafile = os.path.join(monitor_dir, filename)
                    if not os.path.exists(datafile):
                        datafile = os.path.join(self.data_dir, filename)
                    varname = name_pattern.format(n)
                    if os.path.exists(datafile):
                        info_formatted = copy.deepcopy(info)
                        info_formatted["pattern"] = info_formatted[
                            "pattern"].format(fluid_name, "{}")
                        info_formatted["datafile"] = datafile
                        fieldLoader = load1d.FieldLoader1d(varname, info_formatted,
                                                    self)
                        fl.register_variable(varname, "1d", fieldLoader)
                    if not "for each planet" in info or not info[
                            "for each planet"]:
                        break

    def get_scalars(self):
        for fluid_name in self.fluids:
            fl = self.fluids[fluid_name]
            monitor_dir = os.path.join(self.data_dir, "monitor", fluid_name)
            for name_pattern, info in defs.vars_scalar.items():
                for n in range(len(self.planets)):
                    datafile = os.path.join(monitor_dir,
                                            info["file"].format(n))
                    varname = name_pattern.format(n)
                    if os.path.exists(datafile):
                        fl.register_variable(
                            varname, "scalar",
                            loadscalar.ScalarLoader(varname, datafile, info, self))
                    if not "for each planet" in info or not info[
                            "for each planet"]:
                        break

    def get_domain_size(self):
        self.Nphi, self.Nr, self.Ntheta = loadgrid.loadNcells(self.data_dir)

    def calc_dimension(self):
        """ Determine the dimension of the output data from the number of cells."""
        if self.Nphi > 2 and self.Nr > 2 and self.Ntheta > 2:
            self.dim = 3
        elif self.Nphi > 2 and self.Nr > 2:
            self.dim = 2
        else:
            self.dim = 1

    def load_times(self):
        self.output_times = loadtime.loadCoarseOutputTimes(self.data_dir,
                                                  self.units["time"])
        self.fine_output_times = loadtime.loadFineOutputTimes(self.data_dir,
                                                     self.units["time"])

    def get_output_time(self, n):
        return self.output_times[n]

    def get_fine_output_time(self, n):
        rv = self.fine_output_times[n]
        return rv

    def get_units(self):
        self.units = loadunits.loadUnits(self.data_dir)


