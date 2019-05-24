code_info = ( "fargocpt", "0.1", "testloader")

import os
import re
import numpy as np
import astropy.units as u
from . import interface
from .. import fluid
from .. import field
from .. import grid
from .. import vector
from .. import particles

def identify(path):
    identifiers = ["misc.dat", "fargo", "Quantities.dat"]
    Nid = len(identifiers)
    seen_ids = 0
    for root, dirs, files in os.walk(path):
        seen_ids += len([ s in files for s in identifiers])
        if Nid == seen_ids:
            return True
    return False

vars2d = {
    "dens" : {
        "pattern" : "gasdens{}.dat",
        "unit" : u.g/u.cm**2
        }
    ,"energy" : {
        "pattern" : "gasenergy{}.dat",
        "unit" : u.erg
        }
    ,"vrad" : {
        "pattern" : "gasvrad{}.dat",
        "unit" : u.cm/u.s,
        "interfaces" : ["r"]
        }
    ,"vazimuth" : {
        "pattern" : "gasvtheta{}.dat",
        "unit" : u.cm/u.s,
        "interfaces" : ["phi"]
        }

    }

vars1d = {
    "mass" : ""
    }

def var_in_files(varpattern, files):
    p = re.compile(varpattern.replace(".", "\.").format("\d+"))
    for f in files:
        if re.match(p, f):
            return True
    return False

def load_scalar(file, var):
    return [1,1]

def get_data_dir(path):
    rv = None
    for root, dirs, files in os.walk(path):
        if "misc.dat" in files:
            rv = root
            break
    if rv is None:
        raise FileNotFoundError("Could not find identifier file 'misc.dat' in any subfolder of '{}'".format(path))
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
                col, name, unitstr = [s.strip() for s in line[len(identifier):].split("|")]
                found_variables[name] = (col, unitstr)
    return found_variables

def load_text_data_file(filepath, varname):
    # get data
    variables = load_text_data_variables(filepath)
    col = variables[varname][0]
    unit = u.Unit(variables[varname][1])
    data = np.genfromtxt(filepath, usecols=int(col))*unit
    return data


class Loader(interface.Interface):

    code_info = code_info

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data_dir = get_data_dir(self.path)
        self.output_times = []

    def scout(self):
        self.get_domain_size()
        self.get_units()
        self.get_fluids()
        self.get_nbodysystems()
        self.get_fields()
        self.get_vectors()
        self.get_nbodysystems()
        self.get_nbody_planet_variables()

    def get_nbodysystems(self):
        planet_ids = []
        p = re.compile("bigplanet(\d).dat")
        for s in os.listdir(self.data_dir):
            m = re.match(p, s)
            if m:
                planet_ids.append(m.groups()[0])
        self.particlegroups["planets"] = particles.NbodySystem("planets")
        planets = self.particlegroups["planets"]
        planets.register_particles(sorted(planet_ids))
        self.planets = planets

    def get_nbody_planet_variables(self):
        planet_variables = load_text_data_variables(os.path.join(self.data_dir,"bigplanet1.dat"))
        for varname in planet_variables:
            self.particlegroups["planets"].register_variable( varname, PlanetVectorLoader( varname, {"datafile_pattern" : os.path.join(self.data_dir, "bigplanet{}.dat")}, self) )
        # register position
        if all([s in planet_variables for s in ["x", "y"]]):
            self.particlegroups["planets"].register_variable( "position", PlanetVectorLoader( "position",
                                    {"datafile_pattern" : os.path.join(self.data_dir, "bigplanet{}.dat"),
                                     "axes" : { "x" : "x", "y" : "y" }}, self) )
        # register velocities
        if all([s in planet_variables for s in ["vx", "vy"]]):
            self.particlegroups["planets"].register_variable( "velocity", PlanetVectorLoader( "velocity",
                                    {"datafile_pattern" : os.path.join(self.data_dir, "bigplanet{}.dat"),
                                     "axes" : { "x" : "vx", "y" : "vy" }}, self) )

    def get_fluids(self):
        self.fluids["gas"] = fluid.Fluid("gas")

    def get_fields(self):
        gas = self.fluids["gas"]
        files = os.listdir(self.data_dir)
        for varname, info in vars2d.items():
            if var_in_files(info["pattern"], files):
                gas.register_variable(varname, "field", FieldLoader(varname, info, self))

    def get_vectors(self):
        gas = self.fluids["gas"]
        datafile = os.path.join(self.data_dir, "Quantities.dat")
        variables = load_text_data_variables(datafile)
        for varname, (column, unitstr) in variables.items():
            gas.register_variable(varname, "vector", VectorLoader(varname, {"datafile" : datafile }, self))
        # add multi axis vectors
        if all([v in variables for v in ["radial kinetic energy", "azimuthal kinetic energy"]]):
            gas.register_variable("kinetic energy", "vector", VectorLoader(varname, {"datafile" : datafile, "axes" : { "r" : "radial kinetic energy", "phi" : "azimuthal kinetic energy" } }, self))

    def get_domain_size(self):
        self.Nr, self.Nphi = np.genfromtxt(os.path.join(self.data_dir, "dimensions.dat"), usecols=(4,5), dtype=int)

    def get_output_time(self, n):
        if self.output_times == []:
            self.output_times = load_text_data_file(os.path.join(self.data_dir, "misc.dat"), "physical time")
        return self.output_times[n]

    def get_units(self):
        with open(os.path.join(self.data_dir,'units.dat'),'r') as f:
            self.units = {l[0] : float(l[1])*u.Unit(l[2]) for l in
                          [l.split() for l in f if l.split()[0] != '#' and len(l.split())==3]}


class FieldLoader:

    def __init__(self, name, info, loader, *args, **kwargs):
        self.loader = loader
        self.info = info
        self.name = name


    def __call__(self, n):
        t = self.loader.get_output_time(n)
        f = field.Field(self.load_grid(n), self.load_data(n), t, self.name)
        return f

    def load_data(self, n):
        if "unit" in self.info:
            unit = self.info["unit"]
        else:
            unit = 1
            for baseunit, power in self.info["unitpowers"].items():
                unit = unit*units[baseunit]**power
        Nr = self.loader.Nr + (1 if "interfaces" in self.info and "r" in self.info["interfaces"] else 0)
        Nphi = self.loader.Nphi + (1 if "interfaces" in self.info and "phi" in self.info["interfaces"] else 0)

        rv = np.fromfile(self.loader.data_dir + "/" + self.info["pattern"].format(n)).reshape(Nr, Nphi)*unit
        return rv

    def load_grid(self, n):
        r_i = np.genfromtxt(self.loader.data_dir + "/used_rad.dat")*self.loader.units["length"]
        phi_i = np.linspace(-np.pi, np.pi, self.loader.Nphi+1)*u.Unit(1)
        active_interfaces = self.info["interfaces"] if "interfaces" in self.info else []
        g = grid.PolarGrid(r_i = r_i, phi_i = phi_i, active_interfaces=active_interfaces)
        return g

class VectorLoader:

    def __init__(self, name, info, loader, *args, **kwargs):
        self.loader = loader
        self.info = info
        self.name = name

    def __call__(self):
        axes = [] if not "axes" in self.info else [key for key in self.info["axes"]]
        f = vector.Vector(self.load_time(), self.load_data(), name=self.name, axes=axes)
        return f

    def load_data(self):
        if "axes" in self.info:
            rv = []
            unit = None
            for ax, varname in self.info["axes"].items():
                d = load_text_data_file(self.info["datafile"], varname)
                if unit is None:
                    unit = d.unit
                else:
                    if unit != d.unit:
                        raise ValueError("Units of multiaxis vector don't match")
                rv.append(d)
            rv = np.array(rv).transpose()*unit
        else:
            rv = load_text_data_file(self.info["datafile"], self.name)
        return rv

    def load_time(self):
        rv = load_text_data_file(self.info["datafile"], "physical time")
        return rv

class PlanetVectorLoader(VectorLoader):

    def __call__(self, num_output, particle_ids):
        axes = [] if not "axes" in self.info else [key for key in self.info["axes"]]
        rv = particles.ParticleVector(self.load_time(particle_ids), self.load_data_multiple(particle_ids), name=self.name, axes=axes)
        return rv

    def load_time(self, particle_ids):
        rv = load_text_data_file(self.info["datafile_pattern"].format(particle_ids[0]), "physical time")
        return rv

    def load_data_multiple(self, particle_ids):
        data = []
        unit = None
        for pid in particle_ids:
            self.info["datafile"] = self.info["datafile_pattern"].format(pid)
            d = self.load_data()
            if unit is None:
                unit = d.unit
            data.append(d.value)

        data = np.array(data)
        data = data*unit
        data = np.moveaxis(data, 0, 1)
        return data
