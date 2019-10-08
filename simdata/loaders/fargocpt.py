code_info = ( "fargocpt", "0.1", "testloader")

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
        seen_ids += len([ 1 for s in identifiers if s in files])
    if seen_ids >=2:
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

alias_fields = {
    "mass density"          : "dens"
    ,"velocity radial"      : "vrad"
    ,"velocity azimuthal"   : "vazimuth"
    ,"total energy density" : "energy"
}

alias_reduced = {
    "output time step"           : "analysis time step"
    ,"simulation time"           : "physical time"
    ,"mass"                      : "mass"
    ,"angular momentum"          : "angular momentum"
    ,"total energy"              : "total energy"
    ,"internal energy"           : "internal energy"
    ,"kinetic energy"            : "kinematic energy"
    ,"kinetic energy radial"     : "radial kinetic energy"
    ,"kinetic energy azimuthal"  : "azimuthal kinetic energy"
    ,"eccentricity"              : "eccentricity"
    ,"argument of periapsis"     : "periastron"
    ,"mass flow inner"           : ""
    ,"mass flow outer"           : ""
    ,"mass flow wavedamping"     : ""
    ,"mass flow densityfloor"    : ""
}

alias_particle = {
    "output time step"  : "time step"
    ,"simulation time"  : "physical time"
    ,"position"         : "position"
    ,"velocity"         : "velocity"
    ,"mass"             : "mass"
    ,"angular momentum" : "angular momentum"
    ,"eccentricity"     : "eccentricity"
    ,"semi-major axis"  : "semi-major axis"
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
        self.get_scalars()
        self.get_planets()
        self.get_nbodysystems()
        self.register_alias()

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
            planet_variables = load_text_data_variables(os.path.join(self.data_dir,"bigplanet{}.dat".format(pid)))
            for varname in planet_variables:
                datafile =  os.path.join(self.data_dir, "bigplanet{}.dat".format(pid))
                loader = ScalarLoader( varname, datafile, self)
                planet.register_variable( varname, loader)


    def get_fluids(self):
        self.fluids["gas"] = fluid.Fluid("gas")

    def get_fields(self):
        self.get_fields_2d()

    def get_fields_2d(self):
        gas = self.fluids["gas"]
        files = os.listdir(self.data_dir)
        for varname, info in vars2d.items():
            if var_in_files(info["pattern"], files):
                gas.register_variable(varname, "2d", FieldLoader2d(varname, info, self))

    def get_scalars(self):
        gas = self.fluids["gas"]
        datafile = os.path.join(self.data_dir, "Quantities.dat")
        variables = load_text_data_variables(datafile)
        for varname, (column, unitstr) in variables.items():
            gas.register_variable(varname, "scalar", ScalarLoader(varname, datafile , self))

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

class FieldLoader2d(interface.FieldLoader):

    def load_time(self, n):
        rv = self.loader.get_output_time(n)
        return rv

    def load_data(self, n):
        if "unit" in self.info:
            unit = self.info["unit"]
        else:
            unit = 1
            for baseunit, power in self.info["unitpowers"].items():
                unit = unit*units[baseunit]**power
        Nr = self.loader.Nr + (1 if "interfaces" in self.info and "r" in self.info["interfaces"] else 0)
        Nphi = self.loader.Nphi #+ (1 if "interfaces" in self.info and "phi" in self.info["interfaces"] else 0)

        rv = np.fromfile(self.loader.data_dir + "/" + self.info["pattern"].format(n)).reshape(Nr, Nphi)*unit
        return rv

    def load_grid(self, n):
        r_i = np.genfromtxt(self.loader.data_dir + "/used_rad.dat")*self.loader.units["length"]
        phi_i = np.linspace(-np.pi, np.pi, self.loader.Nphi+1)*u.Unit(1)
        active_interfaces = self.info["interfaces"] if "interfaces" in self.info else []
        g = grid.PolarGrid(r_i = r_i, phi_i = phi_i, active_interfaces=active_interfaces)
        return g

class ScalarLoader:
    def __init__(self, name, datafile, loader, *args, **kwargs):
        self.loader = loader
        self.datafile = datafile
        self.name = name

    def __call__(self):
        time = self.load_time()
        data = self.load_data()
        f = scalar.Scalar(time, data, name=self.name )
        return f

    def load_data(self):
        rv = load_text_data_file(self.datafile, self.name)
        return rv

    def load_time(self):
        rv = load_text_data_file(self.datafile, "physical time")
        return rv
