import os
import re

import astropy.units as u
import numpy as np
from simdata import fluid, particles
from simdata.loaders import interface

from . import defs, load1d, load2d, loadparams, loadscalar

code_info = ("fargocpt", "0.1", "legacy_output")


def identify(path):
    identifiers = ["misc.dat", "fargo", "Quantities.dat"]
    seen_ids = 0
    for _, _, files in os.walk(path):
        seen_ids += len([1 for s in identifiers if s in files])
        if seen_ids >= 2:
            return True
    return False


def var_in_files(varpattern, files):
    p = re.compile(varpattern.replace(".", r"\.").format(r"\d+"))
    for f in files:
        if re.match(p, f):
            return True
    return False


def get_data_dir(path):
    rv = None
    # guess first
    for guess in ["outputs", "output", "out"]:
        guess_dir = os.path.join(path, guess)
        if os.path.isfile(os.path.join(guess_dir, "misc.dat")):
            rv = guess_dir
            break
    # now search whole dir tree
    if rv is None:
        for root, _, files in os.walk(path):
            if "misc.dat" in files:
                rv = root
                break
    if rv is None:
        raise FileNotFoundError(
            "Could not find identifier file 'misc.dat' in any subfolder of '{}'"
            .format(path))
    return rv


class Loader(interface.Interface):

    code_info = code_info

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data_dir = get_data_dir(self.path)
        self.output_times = []
        self.fine_output_times = []

    def scout(self):
        self.get_units()
        self.get_domain_size()
        self.get_parameters()
        self.load_grid()
        self.load_times()
        self.get_fluids()
        self.get_planets()
        self.get_nbodysystems()
        self.get_fields()
        self.get_scalars()
        self.get_nbodysystems()
        self.register_alias()

    def get_parameters(self):
        param_file = os.path.join(self.data_dir, "../setup/in.par")
        self.parameters = loadparams.get_parameters(param_file)

    def get_domain_size(self):
        self.Nr, self.Nphi = np.genfromtxt(os.path.join(
            self.data_dir, "dimensions.dat"),
            usecols=(4, 5),
            dtype=int)

    def get_units(self):
        with open(os.path.join(self.data_dir, 'units.dat'), 'r') as f:
            self.units = {
                l[0]: float(l[1]) * u.Unit(l[2])
                for l in [
                    l.split() for l in f
                    if l.split()[0] != '#' and len(l.split()) == 3
                ]
            }

    def load_grid(self):
        self.r_i = np.genfromtxt(self.data_dir +
                                 "/used_rad.dat") * self.units["length"]
        self.phi_i = np.linspace(0, 2 * np.pi, self.Nphi + 1) * u.rad

    def load_times(self):
        self.output_times = loadscalar.load_text_data_file(
            os.path.join(self.data_dir, "misc.dat"), "physical time")
        self.fine_output_times = loadscalar.load_text_data_file(
            os.path.join(self.data_dir, "Quantities.dat"), "physical time")

    def get_output_time(self, n):
        return self.output_times[n]

    def get_fine_output_time(self, n):
        rv = self.fine_output_times[n]
        return rv

    def register_alias(self):
        for planet in self.planets:
            planet.alias.register_dict(defs.alias_particle)
        self.fluids["gas"].alias.register_dict(defs.alias_fields)
        self.fluids["gas"].alias.register_dict(defs.alias_reduced)

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
            planet_variables = loadscalar.load_text_data_variables(
                os.path.join(self.data_dir, "bigplanet{}.dat".format(pid)))
            for varname in planet_variables:
                datafile = os.path.join(self.data_dir,
                                        "bigplanet{}.dat".format(pid))
                loader = loadscalar.ScalarLoader(varname, datafile, self)
                planet.register_variable(varname, loader)

    def get_fluids(self):
        self.fluids["gas"] = fluid.Fluid("gas")

    def get_fields(self):
        self.get_fields_2d()
        self.get_fields_1d()

    def get_fields_2d(self):
        gas = self.fluids["gas"]
        files = os.listdir(self.data_dir)
        for varname, info in defs.vars2d.items():
            if var_in_files(info["pattern"], files):
                loader = load2d.FieldLoader2d(varname, info, self)
                gas.register_variable(varname, "2d", loader)

    def get_fields_1d(self):
        gas = self.fluids["gas"]
        files = os.listdir(self.data_dir)
        for n in range(len(self.planets)):
            if "gas1D_torque_planet{}_0.dat".format(n) in files:
                varname = "torque planet {}".format(n)
                path_pattern = os.path.join(
                    self.data_dir,
                    "gas1D_torque_planet{}_{}.dat".format(n, "{}"))
                info = {
                    "pattern": path_pattern
                }
                loader = load1d.FieldLoader1dTorq(varname, info, self)
                gas.register_variable(varname, "1d", loader)
        if "gasMassFlow1D.info" in files:
            varname = "mass flow"
            info = {}
            loader = load1d.FieldLoader1dMassFlow(varname, info, self)
            gas.register_variable(varname, "1d", loader)
        for fname in files:
            m = re.search(r"(.*)1D\.info", fname)
            if m:
                stem = m.groups()[0]
                infofile = os.path.join(self.data_dir, fname)
                varname = stem[3:].lower()
                if varname == "massflow":
                    continue
                loader = load1d.FieldLoader1d(
                    varname, {"infofile": infofile}, self)
                gas.register_variable(varname, "1d", loader)

    def get_scalars(self):
        gas = self.fluids["gas"]
        datafile = os.path.join(self.data_dir, "Quantities.dat")
        variables = loadscalar.load_text_data_variables(datafile)
        for varname, _ in variables.items():
            loader = loadscalar.ScalarLoader(varname, datafile, self)
            gas.register_variable(varname, "scalar", loader)
