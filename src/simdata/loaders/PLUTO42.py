code_info = ("PLUTO", "4.2", "vanilla")
scalar_filename = "info.dat"
log_file = 'pluto.log'
log_file_backup = '../log.log'

import os
import re
import copy
import subprocess  #needed for a call to "tail" in assign_variables()
import numpy as np
import astropy.units as u
import astropy.constants as const
from . import interface
from .. import fluid
from .. import field
from .. import grid
from .. import scalar
from .. import particles
from .PLUTO42_aux import PLUTOgrid as PGrid


def identify(path):
    try:
        get_data_dir(path)
        return True
    except FileNotFoundError:
        return False


NOT_FOUND = -404
CAN_CONVERT = -1

vars_2d = {
    "mass density": {
        "pattern": "rho.{}.dbl",
        "shorthand": "rho",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "mass": 1,
            "length": -2
        }
    },
    "pressure": {
        "pattern": "prs.{}.dbl",
        "shorthand": "prs",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "mass": 1,
            "time": -2
        }
    },
    "vrad": {
        "pattern": "vx1.{}.dbl",
        "shorthand": "vx1",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "length": 1,
            "time": -1
        }
    },
    "vazimuth": {
        "pattern": "vx3.{}.dbl",
        "shorthand": "vx3",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "length": 1,
            "time": -1
        }
    },
    "vpolar": {
        "pattern": "vx2.{}.dbl",
        "shorthand": "vx2",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "length": 1,
            "time": -1
        }
    },
    "vx": {
        "pattern": "vx1.{}.dbl",
        "shorthand": "vx1",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "length": 1,
            "time": -1
        }
    },
    "vy": {
        "pattern": "vx2.{}.dbl",
        "shorthand": "vx2",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "length": 1,
            "time": -1
        }
    },
    "vz": {
        "pattern": "vx3.{}.dbl",
        "shorthand": "vx3",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "length": 1,
            "time": -1
        }
    },
    "temperature": {
        "pattern": "T.{}.dbl",
        "shorthand": "T",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "temperature": 1
        }
    },
    "vortensity": {
        "pattern": "vort.{}.dbl",
        "shorthand": "vort",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "mass": -1,
            "length": 2,
            "time": -1
        }
    },
    "opacity": {
        "pattern": "opac.{}.dbl",
        "shorthand": "opac",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "mass": -1,
            "length": 2
        }
    },
    "pressure scale height": {
        "pattern": "H.{}.dbl",
        "shorthand": "H",
        "numvar": NOT_FOUND,
        "unitpowers": {
            "length": 1
        }
    },
    "aspect ratio": {
        "pattern": "h.{}.dbl",
        "shorthand": "h",
        "numvar": NOT_FOUND,
        "unitpowers": {}
    }
}

planet_vars_scalar = dict()
vars_1d = dict()
#todo
vars_scalar = {
    'mass': {
        'file': scalar_filename,
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {
            "mass": 1
        }
    },
    'thermal energy density': {
        'file': scalar_filename,
        'datacol': 2,
        'timecol': 0,
        'unitpowers': {
            "mass": 1,
            "length": 2,
            "time": -2
        }
    },
    'torque planet': {
        'file': scalar_filename,
        #'for each planet' : True,
        'datacol': 3,
        'timecol': 0,
        'unitpowers': {
            "mass": 1,
            "length": 2,
            "time": -2
        }
    }
}

alias_fields = {
    "velocity radial": "vrad",
    "velocity azimuthal": "vazimuth",
    "thermal energy density": "energy density"
}

alias_reduced = {
    "output time step": "analysis time step",
    "simulation time": "physical time",
    "mass": "mass",
    "angular momentum": "angular momentum",
    "total energy": "total energy",
    "internal energy": "internal energy",
    "kinetic energy": "kinetic energy",
    "eccentricity": "eccentricity",
    "periastron": "periastron",
    "mass flow inner": "",
    "mass flow outer": "",
    "mass flow wavedamping": "",
    "mass flow densityfloor": ""
}

alias_particle = {
    "output time step": "time step",
    "simulation time": "physical time",
    "argument of periapsis": "argument of periapsis",
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
    ptrn = re.compile("grid.out")
    for root, dirs, files in os.walk(path):
        for f in files:
            m = re.search(ptrn, f)
            if m:
                rv = root
                break
    if rv is None:
        raise FileNotFoundError(
            "Could not find identifier file 'grid.out' in any subfolder of '{}'"
            .format(path))
    return rv


def get_unit_from_powers(unitpowers, units):
    unit = 1.0
    for u, p in unitpowers.items():
        unit = unit * units[u]**p
    return unit


class Loader(interface.Interface):

    code_info = code_info

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data_dir = get_data_dir(self.path)
        self.output_times = np.array([])
        self.fine_output_times = np.array([])

    def set_data_dir(self, data_dir):
        self.data_dir = data_dir

    def scout(self):
        self.get_geometry()
        self.get_domain_size()
        self.get_units()
        self.get_fluids()
        self.assign_variables()
        self.assign_convertible_variables()
        self.set_EOS()
        self.set_parameters()
        self.apply_units()
        self.load_times()
        self.get_fields()
        self.get_scalars()
        self.register_alias()

    def apply_units(self):
        for vardict in [planet_vars_scalar, vars_2d, vars_1d, vars_scalar]:
            for var, info in vardict.items():
                info["unit"] = get_unit_from_powers(info["unitpowers"],
                                                    self.units)

    def register_alias(self):
        for particlegroup in self.particlegroups:
            particlegroup.alias.register_dict(alias_particle)
        for planet in self.planets:
            planet.alias.register_dict(alias_particle)
        for name, fluid in self.fluids.items():
            fluid.alias.register_dict(alias_fields)
            fluid.alias.register_dict(alias_reduced)

    def get_fluids(self):
        self.fluids['gas'] = fluid.Fluid('gas')  # hack

    def set_parameters(self, mean_mol_weight=2.353, mass_sys=1.0):
        self.fluids['gas'].mean_mol_weight = mean_mol_weight * u.g / u.mol
        self.mass_sys = mass_sys * self.units['mass'].cgs
        self.GM = (mass_sys * self.units['mass'] * const.G).decompose().cgs
        self.fluids['gas'].mu_R = (self.fluids['gas'].mean_mol_weight /
                                   const.R).decompose().cgs

    def set_EOS(self):
        if vars_2d['pressure']["numvar"] >= 0:
            self.eos = 'IDEAL'
        else:
            self.eos = 'ISOTHERMAL'

    def get_parameters(self):
        mu = self.fluids['gas'].mean_mol_weight
        M_sys = self.mass_sys

        return dict(mean_mol_weight=mu, mass_sys=M_sys)

    def get_fields(self):
        self.get_fields_2d()
        self.get_fields_1d()

    def get_domain_size(self):
        grid = [PGrid.Grid(DIR=i, origin=self.data_dir) for i in [0, 1, 2]]
        self.NX1, self.NX2, self.NX3 = [g.Ncells for g in grid]

    def get_geometry(self):
        self.dimensions, self.geometry, self.coordinates = PGrid.resolve_geometry(
            origin=self.data_dir)

    def assign_variables(self):
        #warning: the following lines only work with python3+
        #get number of columns of file
        with open(self.data_dir + '/dbl.out', 'r') as dblout:
            line = dblout.readline().strip().split()
        _, _, _, _, output_format, _, *varnames = line

        self.output_format = output_format

        #format of varnames: [rho, vx1, vx2?, vx3?, prs?, *additionalvars]
        #vx2, vx3 only exist if there's enough dimensions
        #prs only exists if the simulation is NOT isothermal

        #by default, I expect either (x; y; z) or (r; theta; phi)
        if self.geometry == 'POLAR':  #aka (r; phi; z) -> move phi to 2nd index
            vars_2d["vazimuth"]["pattern"] = "vx2.{}.dbl"
            vars_2d["vazimuth"]["shorthand"] = "vx2"
        elif self.geometry == 'CYLINDRICAL':  #aka (r; z) -> move z to 2nd index
            vars_2d["vz"]["pattern"] = "vx2.{}.dbl"
            vars_2d["vz"]["shorthand"] = "vx2"

        if self.dimensions == 3:
            vars_2d["mass density"]["unitpowers"]["length"] -= 1
            vars_2d["pressure"]["unitpowers"]["length"] -= 1
            vars_scalar["thermal energy density"]["unitpowers"]["length"] -= 1

        #rho is always first
        vars_2d["mass density"]["numvar"] = 0
        #assign velocity components
        for i, direction in enumerate(self.coordinates):
            vars_2d['v' + direction]["numvar"] = i + 1

        #this index will be useful when checking for additional variables
        last_index = self.dimensions + 1

        #check if the simulation is adiabatic
        if 'prs' in varnames:
            vars_2d['pressure']["numvar"] = last_index
            last_index += 1

        #activate additional variables
        for varname in varnames[last_index:]:
            for var in vars_2d:
                if vars_2d[var]["shorthand"] == varname:
                    vars_2d[var]["numvar"] = last_index
                    last_index += 1

        self.NVAR = last_index

        #check if all fields exist in the same file:
        if output_format == 'single_file':
            for var in vars_2d:
                vars_2d[var]["pattern"] = "data.{:04d}.dbl"

    def assign_convertible_variables(self):
        """This is aimed at a future module of simdata.
        PLUTO only outputs a handful of fields, but the user might want to calculate some additional ones.
        This function checks if the conditions are met to do so.
        E.g. pressure and density are needed for temperature, but temperature is not always output."""
        have_prs = vars_2d["pressure"]["numvar"] >= 0
        have_rho = vars_2d["mass density"]["numvar"] >= 0
        have_tmp = False

        if vars_2d["temperature"][
                "numvar"] == NOT_FOUND and have_prs and have_rho:
            vars_2d["temperature"]["numvar"] = CAN_CONVERT
            have_tmp = True

        if have_rho and have_tmp and vars_2d["opacity"]["numvar"] == NOT_FOUND:
            vars_2d["opacity"]["numvar"] = CAN_CONVERT

        if have_tmp:
            vars_2d["pressure scale height"]["numvar"] = CAN_CONVERT
            vars_2d["aspect ratio"]["numvar"] = CAN_CONVERT

    def get_fields_2d(self):
        for fluidname in self.fluids.keys():
            for varname, info in vars_2d.items():
                if vars_2d[varname]["numvar"] != NOT_FOUND:
                    fieldLoader = FieldLoader2d(varname, info, self)
                    self.fluids[fluidname].register_variable(
                        varname, '2d', fieldLoader)

    def get_fields_1d(self):
        pass

    def get_scalars(self):
        #set scalar_filename at top of file!
        for fluid_name in self.fluids:
            fl = self.fluids[fluid_name]
            datafile = os.path.join(self.data_dir, scalar_filename)
            if os.path.exists(datafile):
                for varname, info in vars_scalar.items():
                    fl.register_variable(
                        varname, "scalar",
                        ScalarLoader(varname, datafile, info, self))

    def load_times(self):
        self.output_times = loadCoarseOutputTimes(self.data_dir,
                                                  self.units["time"])
        self.fine_output_times = loadFineOutputTimes(self.data_dir,
                                                     self.units["time"])

    def get_output_time(self, n):
        return self.output_times[n]

    def get_fine_output_time(self, n):
        rv = self.fine_output_times[n]
        return rv

    def get_units(self):
        self.units = loadUnits(self.data_dir, dimensions=self.dimensions)


class FieldLoader2d(interface.FieldLoader):
    def load_time(self, n):
        rv = self.loader.get_output_time(n)
        return rv

    def load_data(self, n):
        """
        simdata expects output in the form data[x1][x2][x3], but most operations (for disks!)
        are much easier if the x1 coordinate (r) is exposed first (i.e. data[x3][x2][x1]).
        This function returns data in the [x1..3] format, but internally it uses [x3..1]
        since I want to call this recursively sometimes.
        This is controlled by setting internal to True for internal use.
        """
        if self.info["numvar"] >= 0:  #a regular case
            rv = self.read_data(n)
        elif self.info["numvar"] == CAN_CONVERT:
            rv = self.convert_data(n)
        else:
            raise KeyError("I could not identify this variable!")

        rv = rv.swapaxes(0, -1)

        return rv

    def load_grid(self, n=0):
        g = [PGrid.Grid(i, origin=self.loader.data_dir) for i in [0, 1, 2]]

        r_i = g[0].xi * self.loader.units["length"]
        phi_i = g[1].xi * u.radian
        return grid.PolarGrid(r_i=r_i, phi_i=phi_i, active_interfaces=[])

    def read_data(self, n, qty=None):
        """Given a file number n, reads the corresponding data file and returns the data.

        This function is meant to be used internally only!
        It's the backend for load_data().

        The argument "qty" is used CURRENTLY INTERNALLY ONLY.
        In the future there will be a conversion module that might depend on it.
        By calling with qty!=None, the self.info argument is overridden
        and the requested quantity is returned itself.

        Issue: qty is searched for in vars_2d, not registered variables.
        """
        file_format = [self.loader.NX3, self.loader.NX2,
                       self.loader.NX1][3 - self.loader.dimensions:]

        if not qty: qty_info = self.info
        else: qty_info = vars_2d[qty]

        unit = qty_info["unit"]

        filename = self.loader.data_dir + "/" + qty_info["pattern"].format(n)
        if self.loader.output_format == 'single_file':
            memmap = np.memmap(filename,
                               dtype=float,
                               shape=(self.loader.NVAR, *file_format))
            rv = np.array(memmap[qty_info["numvar"]], dtype=float)
            rv = (rv * unit).decompose().cgs
        else:
            rv = np.fromfile(filename).reshape(*file_format)
            rv = (rv * unit).decompose().cgs

        return rv

    def convert_data(self, n, qty=None):
        """Given a requested field (e.g. "temperature"), calculates said field using existing data.
        For temperature, it will read pressure and density and do the math needed.
        """
        if not qty: name = self.name
        else: name = qty

        if name == 'temperature':
            mu_R = self.loader.fluids['gas'].mu_R
            prs = self.read_data(n, 'pressure')
            dens = self.read_data(n, 'mass density')
            rv = (prs / dens * mu_R).decompose().cgs

        elif name == 'pressure scale height':
            prs = self.read_data(n, 'pressure')
            dens = self.read_data(n, 'mass density')
            GM = self.loader.GM
            r = self.load_grid().r_c
            rv = np.sqrt(prs / (dens * GM) * (r * r * r)).decompose().cgs

        elif name == 'aspect ratio':
            H = self.convert_data(n, 'pressure scale height')
            r = self.load_grid().r_c
            rv = (H / r).decompose().cgs

        return rv


class ScalarLoader:
    def __init__(self, name, datafile, info, loader, *args, **kwargs):
        self.loader = loader
        self.datafile = datafile
        self.info = info
        self.name = name
        self.units = loader.units

    def __call__(self):
        time = self.load_time()
        data = self.load_data()
        f = scalar.Scalar(time, data, name=self.name)
        return f

    def load_data(self):
        col = self.info["datacol"]
        unit = self.info["unit"]
        rv = np.genfromtxt(self.datafile, usecols=int(col),
                           skip_header=1) * unit
        return rv

    def load_time(self):
        col = self.info["timecol"]
        unit = self.units["time"]
        rv = np.genfromtxt(self.datafile, usecols=int(col),
                           skip_header=1) * unit
        return rv


def loadCoarseOutputTimes(dataDir, unit):
    timestamps = np.genfromtxt(dataDir + '/dbl.out',
                               usecols=(1),
                               unpack=True,
                               dtype=float)

    return timestamps * unit


def loadFineOutputTimes(dataDir, unit):
    timestamps = np.genfromtxt(dataDir + '/' + scalar_filename,
                               usecols=(0),
                               unpack=True,
                               skip_header=1,
                               dtype=float)

    return timestamps * unit


#todo, call grid
def loadRadius(dataDir, unit, interface=False):
    g = PGrid.Grid(DIR=0, origin=self.loader.data_dir)

    if not interface: return (g.x * unit, g.dx * unit)
    else: return (g.xi * unit, g.dx * unit)


#todo, call grid
def loadPhi(dataDir, interface=False):
    g = PGrid.Grid(DIR=0, origin=self.loader.data_dir)

    if not interface: return g.x
    else: return g.xi


def loadMeshGrid(dataDir, unit):
    # return a meshgrid for the disk to plot data
    R, Phi = loadMeshGridPolar(dataDir, unit)
    X = R * np.cos(Phi)
    Y = R * np.sin(Phi)
    return (X, Y)


def loadMeshGridPolar(dataDir, unit):
    phi = loadPhi(dataDir)
    r, dr = loadRadius(dataDir, unit)
    Phi, R = np.meshgrid(phi, r)
    return (R, Phi)


def loadUnits(dataDir, dimensions):
    # extract lines for units from pluto log file (specify at top of file!)
    unit_lines = []
    in_unit_block = False
    if os.path.isfile(os.path.join(dataDir, log_file)): use_log_file = log_file
    else: use_log_file = log_file_backup

    with open(os.path.join(dataDir, use_log_file)) as f:
        for line in f:
            line = line.strip()
            if line == "> Normalization Units:":
                in_unit_block = True
                continue
            if in_unit_block:
                if line == "":
                    continue
                elif line[0] == ">":
                    break
                else:
                    unit_lines.append(line)
    # parse lines into units dict
    units = {}
    for line in unit_lines:
        parts = line.split()
        name = parts[0].strip(":[]").lower()
        value = float(parts[1])
        unit_str = parts[2].strip(",()")
        if unit_str == "X":
            unit_str = parts[-1].strip("()")
        unit = u.Unit(unit_str.replace('gr', 'g').replace('sec', 's'))
        units[name] = u.Unit(value * unit)
    # calculate mass unit from density and length
    units["density"] *= units["length"]**(3 - dimensions
                                          )  # fix for non-3d density
    units["mass"] = u.Unit(
        (units["density"] * units["length"]**dimensions).decompose()).cgs
    units["time"] = units["time"] / (1 * units["length"]).value  #why?
    return units


def loadNcells(dataDir):
    #todo
    # Nphi, Nr = np.genfromtxt(os.path.join(dataDir, 'dimensions.dat'), usecols=(6,7), dtype=int)
    grid = [PGrid.Grid(DIR=i, origin=dataDir) for i in [0, 1, 2]]
    NX1, NX2, NX3 = [g.Ncells for g in grid]

    return (NX1, NX2, NX3)
