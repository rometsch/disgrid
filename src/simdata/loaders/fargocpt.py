code_info = ("fargocpt", "0.1", "testloader")

import os
import re
import numpy as np
import astropy.units as u
from . import interface
from .. import fluid
from .. import field
from .. import grid
from .. import scalar
from .. import particles


def identify(path):
    identifiers = ["misc.dat", "fargo", "Quantities.dat"]
    Nid = len(identifiers)
    seen_ids = 0
    for root, dirs, files in os.walk(path):
        seen_ids += len([1 for s in identifiers if s in files])
    if seen_ids >= 2:
        return True
    return False


vars2d = {
    "dens": {
        "pattern": "gasdens{}.dat",
        "unit": u.g / u.cm**2
    },
    "energy": {
        "pattern": "gasenergy{}.dat",
        "unit": u.erg
    },
    "vrad": {
        "pattern": "gasvrad{}.dat",
        "unit": u.cm / u.s,
        "interfaces": ["r"]
    },
    "vazimuth": {
        "pattern": "gasvtheta{}.dat",
        "unit": u.cm / u.s,
        "interfaces": ["phi"]
    }
}

alias_fields = {
    "mass density": "dens",
    "velocity radial": "vrad",
    "velocity azimuthal": "vazimuth",
    "total energy density": "energy"
}

alias_reduced = {
    "output time step": "analysis time step",
    "simulation time": "physical time",
    "mass": "mass",
    "angular momentum": "angular momentum",
    "total energy": "total energy",
    "internal energy": "internal energy",
    "kinetic energy": "kinematic energy",
    "kinetic energy radial": "radial kinetic energy",
    "kinetic energy azimuthal": "azimuthal kinetic energy",
    "eccentricity": "eccentricity",
    "argument of periapsis": "periastron",
    "mass flow inner": "",
    "mass flow outer": "",
    "mass flow wavedamping": "",
    "mass flow densityfloor": ""
}

alias_particle = {
    "output time step": "time step",
    "simulation time": "physical time",
    "position": "position",
    "velocity": "velocity",
    "mass": "mass",
    "angular momentum": "angular momentum",
    "eccentricity": "eccentricity",
    "semi-major axis": "semi-major axis"
}


def var_in_files(varpattern, files):
    p = re.compile(varpattern.replace(".", "\.").format("\d+"))
    for f in files:
        if re.match(p, f):
            return True
    return False


def load_scalar(file, var):
    return [1, 1]


def get_data_dir(path):
    rv = None
    for root, dirs, files in os.walk(path):
        if "misc.dat" in files:
            rv = root
            break
    if rv is None:
        raise FileNotFoundError(
            "Could not find identifier file 'misc.dat' in any subfolder of '{}'"
            .format(path))
    return rv


def load_text_data_variables(filepath):
    # load all variable definitions from a text file
    # which contains the variable names and colums in its header.
    # each variable is indicated by
    # "#variable: {column number} | {variable name} | {unit}"
    found_variables = {}
    with open(os.path.join(filepath)) as f:
        for line in f:
            line = line.strip()
            if line[0] != "#":
                break
            identifier = "#variable:"
            if line[:len(identifier)] == identifier:
                col, name, unitstr = [
                    s.strip() for s in line[len(identifier):].split("|")
                ]
                found_variables[name] = (col, unitstr)
    return found_variables


def load_text_data_file(filepath, varname):
    # get data
    variables = load_text_data_variables(filepath)
    col = variables[varname][0]
    unit = u.Unit(variables[varname][1])
    data = np.genfromtxt(filepath, usecols=int(col)) * unit
    return data


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
        self.load_grid()
        self.load_times()
        self.get_fluids()
        self.get_planets()
        self.get_nbodysystems()
        self.get_fields()
        self.get_scalars()
        self.get_nbodysystems()
        self.register_alias()

    def load_grid(self):
        self.r_i = np.genfromtxt(self.data_dir +
                                 "/used_rad.dat") * self.units["length"]
        self.phi_i = np.linspace(-np.pi, np.pi, self.Nphi + 1) * u.Unit(1)

    def load_times(self):
        self.output_times = load_text_data_file(
            os.path.join(self.data_dir, "misc.dat"), "physical time")
        self.fine_output_times = load_text_data_file(
            os.path.join(self.data_dir, "Quantities.dat"), "physical time")

    def get_output_time(self, n):
        return self.output_times[n]

    def get_fine_output_time(self, n):
        rv = self.fine_output_times[n]
        return rv

    def register_alias(self):
        for planet in self.planets:
            planet.alias.register_dict(alias_particle)
        self.fluids["gas"].alias.register_dict(alias_fields)
        self.fluids["gas"].alias.register_dict(alias_reduced)

    def get_nbodysystems(self):
        pass

    def get_planets(self):
        planet_ids = []
        p = re.compile("bigplanet(\d).dat")
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
            planet_variables = load_text_data_variables(
                os.path.join(self.data_dir, "bigplanet{}.dat".format(pid)))
            for varname in planet_variables:
                datafile = os.path.join(self.data_dir,
                                        "bigplanet{}.dat".format(pid))
                loader = ScalarLoader(varname, datafile, self)
                planet.register_variable(varname, loader)

    def get_fluids(self):
        self.fluids["gas"] = fluid.Fluid("gas")

    def get_fields(self):
        self.get_fields_2d()
        self.get_fields_1d()

    def get_fields_2d(self):
        gas = self.fluids["gas"]
        files = os.listdir(self.data_dir)
        for varname, info in vars2d.items():
            if var_in_files(info["pattern"], files):
                gas.register_variable(varname, "2d",
                                      FieldLoader2d(varname, info, self))

    def get_fields_1d(self):
        gas = self.fluids["gas"]
        files = os.listdir(self.data_dir)
        for n in range(len(self.planets)):
            if "gas1D_torque_planet{}_0.dat".format(n) in files:
                varname = "torque planet {}".format(n)
                info = {
                    "pattern":
                    os.path.join(
                        self.data_dir,
                        "gas1D_torque_planet{}_{}.dat".format(n, "{}"))
                }
                loader = FieldLoader1dTorq(varname, info, self)
                gas.register_variable(varname, "1d", loader)
        if "gasMassFlow1D.info" in files:
            varname = "mass flow".format(n)
            info = {}
            loader = FieldLoader1dMassFlow(varname, info, self)
            gas.register_variable(varname, "1d", loader)
        for fname in files:
            m = re.search("(.*)1D\.info", fname)
            if m:
                stem = m.groups()[0]
                infofile = os.path.join(self.data_dir, fname)
                varname = stem[3:].lower()
                if varname == "massflow":
                    continue
                loader = FieldLoader1d(varname, {"infofile": infofile}, self)
                gas.register_variable(varname, "1d", loader)

    def get_scalars(self):
        gas = self.fluids["gas"]
        datafile = os.path.join(self.data_dir, "Quantities.dat")
        variables = load_text_data_variables(datafile)
        for varname, (column, unitstr) in variables.items():
            gas.register_variable(varname, "scalar",
                                  ScalarLoader(varname, datafile, self))

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


class FieldLoader2d(interface.FieldLoader):
    def load_time(self, n):
        if n is None:
            rv = self.loader.output_times
        else:
            rv = self.loader.get_output_time(n)
        return rv

    def load_data(self, n):
        if "unit" in self.info:
            unit = self.info["unit"]
        else:
            unit = 1
            for baseunit, power in self.info["unitpowers"].items():
                unit = unit * units[baseunit]**power
        Nr = self.loader.Nr + (1 if "interfaces" in self.info
                               and "r" in self.info["interfaces"] else 0)
        Nphi = self.loader.Nphi  #+ (1 if "interfaces" in self.info and "phi" in self.info["interfaces"] else 0)
        rv = np.fromfile(self.loader.data_dir +
                         "/" + self.info["pattern"].format(n)).reshape(
                             Nr, Nphi) * unit
        return rv

    def load_grid(self, n):
        active_interfaces = self.info[
            "interfaces"] if "interfaces" in self.info else []
        g = grid.PolarGrid(r_i=self.loader.r_i,
                           phi_i=self.loader.phi_i,
                           active_interfaces=active_interfaces)
        return g


def parse_1d_info_file(fname):
    data_dir = os.path.dirname(os.path.abspath(fname))
    stem = os.path.basename(os.path.abspath(fname))[:-7]
    pattern = os.path.join(data_dir, "{}1D{}.dat".format(stem, "{}"))
    with open(fname, "r") as infofile:
        for line in infofile:
            line = line.strip()
            parts = [s.strip() for s in line.split("=")]
            if parts[0] == "Nr":
                Nr = int(parts[1])
            elif parts[0] == "unit":
                unit = u.Unit(parts[1])
    return {"unit": unit, "Nr": Nr, "pattern": pattern}


class FieldLoader1d(interface.FieldLoader):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fileinfo = parse_1d_info_file(self.info["infofile"])

    def load_time(self, n):
        if n is None:
            rv = self.loader.output_times
        else:
            rv = self.loader.get_output_time(n)
        return rv

    def load_data(self, n):
        unit = self.fileinfo["unit"]
        rv = load1dRadial(n, self.fileinfo["pattern"], self.fileinfo["unit"])
        return rv

    def load_grid(self, n):
        if self.fileinfo["Nr"] == len(self.loader.r_i):
            active_interfaces = ["r"]
        else:
            active_interfaces = []
        g = grid.PolarGrid(r_i=self.loader.r_i,
                           active_interfaces=active_interfaces)
        return g


class FieldLoader1dTorq(interface.FieldLoader):
    def load_time(self, n):
        rv = self.loader.get_output_time(n)
        return rv

    def load_data(self, n):
        datafile = self.info["pattern"].format(n)
        data = np.fromfile(datafile, dtype=float)
        v = data[1::2] * u.Unit("cm^2 g / s^2")
        return data

    def load_grid(self, n):
        g = grid.PolarGrid(r_i=self.loader.r_i)
        return g


def load1dMassFlow(n, dataDir):
    # parse datafile info
    Nr = 0
    unit = 1
    with open(os.path.join(dataDir, "gasMassFlow1D.info"), "r") as infofile:
        for line in infofile:
            line = line.strip()
            parts = [s.strip() for s in line.split("=")]
            if parts[0] == "Nr":
                Nr = int(parts[1])
            elif parts[0] == "unit":
                unit = u.Unit(parts[1])
    # get data
    Nsamples = 100
    data = np.fromfile(os.path.join(dataDir, "gasMassFlow1D.dat"),
                       count=Nsamples * Nr,
                       offset=Nr * n * 8)
    data = data.reshape([int(len(data) / Nr), Nr])
    data = np.average(data, axis=0)
    if n == 0:
        data = np.zeros(len(data))
    rv = data * unit
    return rv


class FieldLoader1dMassFlow(interface.FieldLoader):
    def load_time(self, n):
        rv = self.loader.get_fine_output_time(n)
        return rv

    def load_data(self, n):
        rv = load1dMassFlow(n, self.data_dir)
        return rv

    def load_grid(self, n):
        r_i = np.genfromtxt(self.loader.data_dir +
                            "/used_rad.dat") * self.loader.units["length"]
        g = grid.PolarGrid(r_i=r_i, active_interfaces=["r"])
        return g


class ScalarLoader:
    def __init__(self, name, datafile, loader, *args, **kwargs):
        self.loader = loader
        self.datafile = datafile
        self.name = name

    def __call__(self):
        time = self.load_time()
        data = self.load_data()
        f = scalar.Scalar(time, data, name=self.name)
        return f

    def load_data(self):
        rv = load_text_data_file(self.datafile, self.name)
        return rv

    def load_time(self):
        rv = load_text_data_file(self.datafile, "physical time")
        return rv


def load1dRadial(n, dataFilePattern, unit, lengthunit=None):
    data = np.fromfile(dataFilePattern.format(n), dtype=float)
    if 'torque_planet' in dataFilePattern:
        v = data[1::2] * unit
        r = data[::2]
    else:
        v = data[1::4] * unit
        r = data[::4]
    if lengthunit is None:
        return v
    else:
        r = r * lengthunit
        return (r, v)


def load1dMassFlow(n, dataDir):
    # parse datafile info
    Nr = 0
    unit = 1
    with open(os.path.join(dataDir, "gasMassFlow1D.info"), "r") as infofile:
        for line in infofile:
            line = line.strip()
            parts = [s.strip() for s in line.split("=")]
            if parts[0] == "Nr":
                Nr = int(parts[1])
            elif parts[0] == "unit":
                unit = u.Unit(parts[1])
    # get data
    Nsamples = 100
    data = np.fromfile(os.path.join(dataDir, "gasMassFlow1D.dat"),
                       count=Nsamples * Nr,
                       offset=Nr * n * 8)
    data = data.reshape([int(len(data) / Nr), Nr])
    data = np.average(data, axis=0)
    if n == 0:
        data = np.zeros(len(data))
    rv = data * unit
    return rv
