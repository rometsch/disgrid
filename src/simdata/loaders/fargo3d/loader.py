import copy
import os
import re

import astropy.constants as const
import astropy.units as u
import numpy as np
from simdata import field, fluid, grid, particles
from simdata.loaders import interface

from . import defs
from . import load1d
from . import load2d
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
    for root, dirs, files in os.walk(path):
        for f in files:
            m = re.search(ptrn, f)
            if m:
                rv = root
                break
    if rv is None:
        raise FileNotFoundError(
            r"Could not find identifier file 'summary\d+.dat' in any subfolder of '{}'"
            .format(path))
    return rv


def get_unit_from_powers(unitpowers, units):
    unit = 1.0
    for k, p in unitpowers.items():
        unit = unit * units[k]**p
    return unit


class Loader(interface.Interface):

    code_info = code_info

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data_dir = get_data_dir(self.path)
        self.output_times = np.array([])
        self.fine_output_times = np.array([])

    def scout(self):
        self.get_domain_size()
        self.get_parameters()
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
        for vardict in [defs.planet_vars_scalar, defs.vars_2d, defs.vars_1d, defs.vars_scalar]:
            for info in vardict.values():
                info["unit"] = get_unit_from_powers(info["unitpowers"],
                                                    self.units)

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
        self.get_fields_2d()
        self.get_fields_1d()

    def get_fields_2d(self):
        files = os.listdir(self.data_dir)
        for fluidname in self.fluids.keys():
            for varname, info in defs.vars_2d.items():
                info_formatted = copy.deepcopy(info)
                info_formatted["pattern"] = info_formatted["pattern"].format(
                    fluidname, "{}")
                if var_in_files(info_formatted["pattern"], files):
                    fieldLoader = load2d.FieldLoader2d(varname, info_formatted, self)
                    self.fluids[fluidname].register_variable(
                        varname, "2d", fieldLoader)

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
        self.Nphi, self.Nr = loadgrid.loadNcells(self.data_dir)

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


