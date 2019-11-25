code_info = ("fargo3d", "2.0", "multifluid")

import os
import re
import copy
import numpy as np
import astropy.units as u
import astropy.constants as const
from . import interface
from .. import fluid
from .. import field
from .. import grid
from .. import scalar
from .. import particles


def identify(path):
    try:
        get_data_dir(path)
        return True
    except FileNotFoundError:
        return False


vars_2d = {
    "mass density": {
        "pattern": "{}dens{}.dat",
        "unitpowers": {
            "mass": 1,
            "length": -2
        }
    },
    "energy density": {
        "pattern": "{}energy{}.dat",
        "unitpowers": {
            "mass": 1,
            "time": -2
        }
    },
    "vrad": {
        "pattern": "{}vy{}.dat",
        "unitpowers": {
            "length": 1,
            "time": -1
        },
        "interfaces": ["r"],
    },
    "vazimuth": {
        "pattern": "{}vx{}.dat",
        "unitpowers": {
            "length": 1,
            "time": -1
        },
        "interfaces": ["phi"],
    },
    "vpolar": {
        "pattern": "{}vz{}.dat",
        "unitpowers": {
            "length": 1,
            "time": -1
        },
        "interfaces": ["theta"],
    },
    "grainsize": {
        "pattern": "{}grainsize{}.dat",
        "unitpowers": {
            "length": 1
        },
    },
    "grainsize drift": {
        "pattern": "{}grainsize_drift{}.dat",
        "unitpowers": {
            "length": 1
        },
    },
    "grainsize frag": {
        "pattern": "{}grainsize_frag{}.dat",
        "unitpowers": {
            "length": 1
        },
    },
    "grainsize driftfrag": {
        "pattern": "{}grainsize_driftfrag{}.dat",
        "unitpowers": {
            "length": 1
        },
    },
    "grainsize coag": {
        "pattern": "{}grainsize_coag{}.dat",
        "unitpowers": {
            "length": 1
        },
    }
}

vars_1d = {
    'torque planet {}': {
        'pattern': 'torq_1d_Y_raw_planet_{}.dat',
        'for each planet': True,
        'directions': ["r"],
        'unitpowers': {
            "mass": 1,
            "length": 2,
            "time": -2
        }
    },
    'velocity radial': {
        'pattern': 'torq_1d_Y_raw_planet_{}.dat',
        'directions': ["r"],
        'unitpowers': {
            "mass": 1,
            "length": 2,
            "time": -2
        }
    },
}

vars_scalar = {
    'mass': {
        'file': 'mass.dat',
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {
            "mass": 1
        }
    },
    'angular momentum': {
        'file': 'momx.dat',
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {
            "mass": 1,
            "length": 2,
            "time": -1
        }
    },
    'kinetic energy azimuthal': {
        'file': 'ekinx.dat',
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {
            "mass": 1,
            "length": 2,
            "time": -2
        }
    },
    'kinetic energy radial': {
        'file': 'ekiny.dat',
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {
            "mass": 1,
            "length": 2,
            "time": -2
        }
    },
    'kinetic energy vertical': {
        'file': 'ekinz.dat',
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {
            "mass": 1,
            "length": 2,
            "time": -2
        }
    },
    'torque planet {}': {
        'file': 'torq_planet_{}.dat',
        'for each planet': True,
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {
            "mass": 1,
            "length": 2,
            "time": -2
        }
    },
}

planet_vars_scalar = {
    'x': {
        'file': 'bigplanet{}.dat',
        'datacol': 1,
        'timecol': 8,
        'unitpowers': {
            'length': 1
        }
    },
    'y': {
        'file': 'bigplanet{}.dat',
        'datacol': 2,
        'timecol': 8,
        'unitpowers': {
            'length': 1
        }
    },
    'z': {
        'file': 'bigplanet{}.dat',
        'datacol': 3,
        'timecol': 8,
        'unitpowers': {
            'length': 1
        }
    },
    'vx': {
        'file': 'bigplanet{}.dat',
        'datacol': 4,
        'timecol': 8,
        'unitpowers': {
            'length': 1,
            "time": -1
        }
    },
    'vy': {
        'file': 'bigplanet{}.dat',
        'datacol': 5,
        'timecol': 8,
        'unitpowers': {
            'length': 1,
            "time": -1
        }
    },
    'vz': {
        'file': 'bigplanet{}.dat',
        'datacol': 6,
        'timecol': 8,
        'unitpowers': {
            'length': 1,
            "time": -1
        }
    },
    'mass': {
        'file': 'bigplanet{}.dat',
        'datacol': 7,
        'timecol': 8,
        'unitpowers': {
            'mass': 1
        }
    },
    'mass': {
        'file': 'bigplanet{}.dat',
        'datacol': 9,
        'timecol': 8,
        'unitpowers': {
            'time': -1
        }
    },
    'time step': {
        'file': 'bigplanet{}.dat',
        'datacol': 0,
        'timecol': 8,
        'unitpowers': {}
    },
    'physical time': {
        'file': 'bigplanet{}.dat',
        'datacol': 8,
        'timecol': 8,
        'unitpowers': {
            "time": 1
        }
    },
    ########################################
    ### orbital elements
    'physical time orbit': {
        'file': 'orbit{}.dat',
        'datacol': 0,
        'timecol': 0,
        'unitpowers': {
            "time": 1
        }
    },
    'eccentricity': {
        'file': 'orbit{}.dat',
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {}
    },
    'semi-major axis': {
        'file': 'orbit{}.dat',
        'datacol': 2,
        'timecol': 0,
        'unitpowers': {
            "length": 1
        }
    },
    'mean anomaly': {
        'file': 'orbit{}.dat',
        'datacol': 3,
        'timecol': 0,
        'unitpowers': {}
    },
    'true anomaly': {
        'file': 'orbit{}.dat',
        'datacol': 4,
        'timecol': 0,
        'unitpowers': {}
    },
    'argument of periapsis': {
        'file': 'orbit{}.dat',
        'datacol': 5,
        'timecol': 0,
        'unitpowers': {}
    },
    'x-axis rotation angle': {
        'file': 'orbit{}.dat',
        'datacol': 6,
        'timecol': 0,
        'unitpowers': {}
    },
    'inclination': {
        'file': 'orbit{}.dat',
        'datacol': 7,
        'timecol': 0,
        'unitpowers': {}
    },
    'ascending node': {
        'file': 'orbit{}.dat',
        'datacol': 8,
        'timecol': 0,
        'unitpowers': {}
    },
    'longitude of periapsis': {
        'file': 'orbit{}.dat',
        'datacol': 9,
        'timecol': 0,
        'unitpowers': {}
    },
}

alias_fields = {
    "velocity radial": "vrad",
    "velocity azimuthal": "vazimuth",
    "total energy density": "energy density"
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
    ptrn = re.compile("summary\d+.dat")
    for root, dirs, files in os.walk(path):
        for f in files:
            m = re.search(ptrn, f)
            if m:
                rv = root
                break
    if rv is None:
        raise FileNotFoundError(
            "Could not find identifier file 'summary\d+.dat' in any subfolder of '{}'"
            .format(path))
    return rv


def find_first_summary(dataDir):
    return "summary{}.dat".format(find_first_summary_number(dataDir))


def find_first_summary_number(dataDir):
    return find_summary_numbers(dataDir)[0]


def find_summary_numbers(dataDir):
    ptrn = re.compile("summary(\d+).dat")
    summaries = []
    for f in os.listdir(dataDir):
        m = re.search(ptrn, f)
        if m:
            n = int(m.groups()[0])
            summaries.append(n)
    summaries.sort()
    return summaries


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
        self.parameters = getParamsFromNthSummary(
            self.data_dir, find_first_summary_number(self.data_dir))

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
            for varname in planet_vars_scalar:
                datafile = os.path.join(self.data_dir,
                                        "bigplanet{}.dat".format(pid))
                loader = ScalarLoader(varname, datafile,
                                      planet_vars_scalar[varname], self)
                planet.register_variable(varname, loader)

    def get_fluids(self):
        ptrn = re.compile("output(.*)\.dat")
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
            for varname, info in vars_2d.items():
                info_formatted = copy.deepcopy(info)
                info_formatted["pattern"] = info_formatted["pattern"].format(
                    fluidname, "{}")
                if var_in_files(info_formatted["pattern"], files):
                    fieldLoader = FieldLoader2d(varname, info_formatted, self)
                    self.fluids[fluidname].register_variable(
                        varname, "2d", fieldLoader)

    def get_fields_1d(self):
        for fluid_name in self.fluids:
            fl = self.fluids[fluid_name]
            monitor_dir = os.path.join(self.data_dir, "monitor", fluid_name)
            monitor_files = os.listdir(monitor_dir)
            for name_pattern, info in vars_1d.items():
                for n in range(len(self.planets)):
                    datafile = os.path.join(monitor_dir,
                                            info["pattern"].format(n))
                    varname = name_pattern.format(n)
                    if os.path.exists(datafile):
                        info_formatted = copy.deepcopy(info)
                        info_formatted["pattern"] = info_formatted[
                            "pattern"].format(fluid_name, "{}")
                        info_formatted["datafile"] = datafile
                        fieldLoader = FieldLoader1d(varname, info_formatted,
                                                    self)
                        fl.register_variable(varname, "1d", fieldLoader)
                    if not "for each planet" in info or not info[
                            "for each planet"]:
                        break

    def get_scalars(self):
        for fluid_name in self.fluids:
            fl = self.fluids[fluid_name]
            monitor_dir = os.path.join(self.data_dir, "monitor", fluid_name)
            monitor_files = os.listdir(monitor_dir)
            for name_pattern, info in vars_scalar.items():
                for n in range(len(self.planets)):
                    datafile = os.path.join(monitor_dir,
                                            info["file"].format(n))
                    varname = name_pattern.format(n)
                    if os.path.exists(datafile):
                        fl.register_variable(
                            varname, "scalar",
                            ScalarLoader(varname, datafile, info, self))
                    if not "for each planet" in info or not info[
                            "for each planet"]:
                        break

    def get_domain_size(self):
        self.Nphi, self.Nr = loadNcells(self.data_dir)

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
        self.units = loadUnits(self.data_dir)


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
        Nr = self.loader.Nr  #+ (1 if "interfaces" in self.info and "r" in self.info["interfaces"] else 0)
        Nphi = self.loader.Nphi  #+ (1 if "interfaces" in self.info and "phi" in self.info["interfaces"] else 0)
        rv = np.fromfile(self.loader.data_dir +
                         "/" + self.info["pattern"].format(n)).reshape(
                             Nr, Nphi) * unit
        return rv

    def load_grid(self, n):
        r_i = np.genfromtxt(self.loader.data_dir + "/domain_y.dat"
                            )[3:-3] * self.loader.units["length"]
        # account for Fargo3d not writing out last radial interface
        if "interfaces" in self.info and "r" in self.info["interfaces"]:
            r_i = r_i[:-1]
        phi_i = np.genfromtxt(self.loader.data_dir +
                              "/domain_x.dat") * u.Unit(1)
        active_interfaces = self.info[
            "interfaces"] if "interfaces" in self.info else []
        g = grid.PolarGrid(r_i=r_i,
                           phi_i=phi_i,
                           active_interfaces=active_interfaces)
        return g


class FieldLoader1d(interface.FieldLoader):
    def load_time(self, n):
        if n is None:
            rv = self.loader.fine_output_times
        else:
            rv = self.loader.get_fine_output_time(n)
        return rv

    def load_data(self, n):
        unit = self.info["unit"]
        Nr = self.loader.Nr  #+ (1 if "interfaces" in self.info and "r" in self.info["interfaces"] else 0)
        Nphi = self.loader.Nphi  #+ (1 if "interfaces" in self.info and "phi" in self.info["interfaces"] else 0)
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
                              "/domain_x.dat") * u.Unit(1)
        active_interfaces = self.info[
            "interfaces"] if "interfaces" in self.info else []
        kwargs = {}
        for d in ["r", "phi"]:
            if d in self.info["directions"]:
                kwargs[d + "_i"] = locals()[d + "_i"]
        kwargs["active_interfaces"] = active_interfaces
        g = grid.PolarGrid(**kwargs)
        return g


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
        rv = np.genfromtxt(self.datafile, usecols=int(col)) * unit
        return rv

    def load_time(self):
        col = self.info["timecol"]
        unit = self.units["time"]
        rv = np.genfromtxt(self.datafile, usecols=int(col)) * unit
        return rv


def loadCoarseOutputTimes(dataDir, unit):
    # search all summary.dat files for the time
    outputTimes = []
    pattern = re.compile('OUTPUT [0-9]* at simulation time ([0-9\.eE+-]*)')
    for f in sorted([f for f in os.listdir(dataDir) if 'summary' in f],
                    key=lambda x: int(x[7:-4])):
        with open(os.path.join(dataDir, f), 'r') as infile:
            datastr = infile.read()
            matches = re.findall(pattern, datastr)
            try:
                outputTimes.append(float(matches[0]))
            except ValueError:
                break
    times = np.array(outputTimes)
    # fall back to reading the planet file for multifluid version
    # which is missing the summary files
    #times = np.genfromtxt( os.path.join(dataDir, 'planet0.dat'))[:,8]
    #times = times*unit
    # correct for double entries in the planet file
    return times * unit


def loadFineOutputTimes(dataDir, unit):
    numbers = find_summary_numbers(dataDir)
    times = np.array([])
    for n in numbers:
        params = getParamsFromNthSummary(dataDir, n)
        dt = params["dt"]
        Ninterm = params["ninterm"]
        offset = 0 if len(times) == 0 else times[-1]
        new_times = np.arange(1, Ninterm + 1) * dt + offset
        times = np.append(times, new_times)
    times = times * unit
    return times


def getParamFromSummary(dataDir, param):
    return getParamsFromNthSummary(
        dataDir, find_first_summary_number(dataDir))[param.lower()]


def getParamsFromNthSummary(dataDir, n):
    # parse the Nth summary file to get all
    search_active = False
    parameters = {}
    with open(os.path.join(dataDir, "summary{}.dat".format(n))) as f:
        for line in f:
            line = line.strip()
            if not search_active:
                # look for the parameter section identifier
                if line == "PARAMETERS SECTION:":
                    search_active = True
                continue
            if line == "" or line[0] in ["#", "="]:
                continue
            if line.startswith("*** Input file: "):
                parameters["config path"] = line.split(":")[-1].strip()
                break
            parts = [s.strip() for s in line.split()]
            try:
                val = int(parts[1])
            except ValueError:
                try:
                    val = float(parts[1])
                except ValueError:
                    val = parts[1]
            parameters[parts[0].lower()] = val
    return parameters


def loadRadius(dataDir, unit, interface=False):
    r = np.genfromtxt(os.path.join(dataDir, 'domain_y.dat')) * unit
    r = r[3:-3]  #remove ghost cells
    dr = r[1:] - r[:-1]
    if not interface:
        r = 0.5 * (r[1:] + r[:-1])
    return (r, dr)


def loadPhi(dataDir, interface=False):
    #phiMin, phiMax, Nphi = np.genfromtxt(os.path.join(dataDir, 'dimensions.dat'), usecols=(0,1,6))
    #phi = np.linspace(phiMin, phiMax, Nphi)
    phi = np.genfromtxt(os.path.join(dataDir, 'domain_x.dat'))
    if not interface:
        phi = 0.5 * (phi[1:] + phi[:-1])
    return phi


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


def loadUnits(dataDir):
    ### load data units
    first_summary = os.path.join(dataDir, find_first_summary(dataDir))
    if os.path.exists(first_summary):
        ptrn = "COMPILATION OPTION SECTION:\n==============================\n.*\-DCGS.*\nGhost"
        with open(first_summary, 'r') as infile:
            if re.search(ptrn, infile.read()):
                # have cgs units
                units = {
                    "mass": u.g,
                    "time": u.s,
                    "length": u.cm,
                    "temperature": u.K
                }
                return units

        # Try to extract unit normalisation from summary
        units = {}
        units["mass"] = 1.0
        units["time"] = 1.0
        units["length"] = 1.0

        with open(first_summary, 'r') as infile:
            ptrn = "(?<=R0 = \()\d+.\d+"
            m = re.search(ptrn, infile.read())
            if m:
                units["length"] *= float(m.group()[0])

        with open(first_summary, 'r') as infile:
            ptrn = "(?<=MSTAR = \()\d+.\d+"
            m = re.search(ptrn, infile.read())
            if m:
                units["mass"] *= float(m.group()[0])

        with open(first_summary, 'r') as infile:
            ptrn = r"STEFANK =.*\*pow\(\(\d+\.\d+\)\/\((\d+\.*\d*\*\d+\.\d*\w+\d*)\),-0\.5"
            m = re.search(ptrn, infile.read())
            if m:
                components = m.group(1).split('*')
                units["length"] *= float(components[0]) * float(
                    components[1]) * u.cm

        with open(first_summary, 'r') as infile:
            ptrn = r"STEFANK =.*\*pow\(\(\d+\.\d*\)\/(\d+\.\d*\w+\d*),-1\.5\)\*"
            m = re.search(ptrn, infile.read())
            if m:
                components = m.group(1).split('*')
                units["mass"] *= float(m.group(1)) * u.g
        units["time"] = (np.sqrt(units["length"]**3 /
                                 (const.G.cgs * units["mass"]))).to(u.s)

    return units
    # now try units file
    # try:
    #     units = {l[0] : float(l[1])*u.Unit(l[2]) for l in
    #              [l.split() for l in open(os.path.join(dataDir,'units.dat'),'r')
    #               if l.split()[0] != '#' and len(l.split())==3]}
    #     ### fix temperature unit
    #     units['temperature'] = 1*u.K
    # except FileNotFoundError:
    # Fall back to default units
    # units = { 'mass' : u.solMass, 'time' : 5.2**1.5*u.yr/(2*np.pi), 'length' : 5.2*u.au }
    # Fall back to dimensionless units
    #units = { bu : 1 for bu in ['mass', 'time', 'length'] }


def loadNcells(dataDir):
    # Nphi, Nr = np.genfromtxt(os.path.join(dataDir, 'dimensions.dat'), usecols=(6,7), dtype=int)
    Nphi = int(getParamFromSummary(dataDir, "Nx"))
    Nr = int(getParamFromSummary(dataDir, "Ny"))
    return (Nphi, Nr)


def load1dRadialMonitorRaw(n, dataFile, Ncells, unit):
    # load data by first seeking the right position
    # and then reading Nrad floats
    f = open(dataFile, "rb")  # reopen the file
    f.seek(n * Ncells * 8, os.SEEK_SET)  # seek
    v = np.fromfile(f, dtype=np.float64, count=Ncells)
    f.close()
    v = v * unit
    return v


def load1dRadialMonitorDensity(n, dataFile, r, dr, unit):
    # Fargo3d outputs monitor variables as the integral over Phi
    # Correct this by computing the density
    rv = load1dRadialMonitorRaw(n, dataFile, len(r), unit)
    rv = rv / np.pi / ((r + dr / 2)**2 - (r - dr / 2)**2)
    return rv


def load1dRadialDensityAveragedFrom2d(n, dataFilePattern, Nr, Nphi, r, dr,
                                      unit):
    rv = load1dRadialAveragedFrom2d(n, dataFilePattern, Nr, Nphi, unit)
    # make it a 2d density
    rv = rv / np.pi / ((r + dr / 2)**2 - (r - dr / 2)**2)
    return rv


def load1dRadialAveragedFrom2d(n, dataFilePattern, Nr, Nphi, unit):
    # Load 2d data and average over the azimuthal domain
    data = load2d(n, dataFilePattern, Nr, Nphi, unit)
    rv = np.mean(data, axis=1)
    return rv


def load2d(n, dataFilePattern, Nr, Nphi, unit):
    # Load 2d data an reshape it
    rv = np.fromfile(dataFilePattern.format(n)).reshape(Nr, Nphi) * unit
    return rv
