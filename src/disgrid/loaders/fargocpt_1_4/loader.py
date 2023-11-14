import os
import re
import yaml

from weakref import proxy

import astropy.units as u
import numpy as np
from ... import fluid, particles
from .. import interface

from . import defs, load1d, load2d, loadparams, loadscalar, loadparticles

identifier = "fargocpt_output_v1_4"
code_info = ("fargocpt", "1.4", "output v1.4")

smoking_gun = identifier

def identify(path):
    
    for _, dirs, files in os.walk(path, topdown=True):
        if "snapshots" in dirs:
            dirs.remove("snapshots")
        if identifier in files:
            return True
    return False


def var_in_files(varpattern, files):
    p = re.compile(varpattern.replace(".", r"\.").format(r"\d+"))
    for f in files:
        if re.match(p, f):
            return True
    return False


def quantity_from_spec(spec):
    """ Generate a astropy quantity object from a (value, unit string) tuple.

    Parameters
    ----------
    spec : (array, str)
        Tuple storing the values and the unit string.
    """
    return np.array(spec[0])*u.Unit(spec[1])


def get_data_dir(path):
    rv = None
    # guess first
    for guess in ["outputs", "output", "out"]:
        guess_dir = os.path.join(path, guess)
        if os.path.isfile(os.path.join(guess_dir, identifier)):
            rv = guess_dir
            break
    # now search whole dir tree
    if rv is None:
        for root, dirs, files in os.walk(path, topdown=True):
            if "snapshots" in dirs:
                dirs.remove("snapshots")
            if identifier in files:
                rv = root
                break
    if rv is None:
        raise FileNotFoundError(
            f"Could not find identifier file '{identifier}' in any subfolder of '{path}'")
    return rv


class Loader(interface.Interface):

    code_info = code_info

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.update:
            self.uptodate = False
        elif "spec" in kwargs and kwargs["spec"] is not None:
            self.spec = kwargs["spec"].copy()
            self.data_dir = os.path.join(self.path, self.spec["data_dir"])
            try:
                self.get_last_activity()
                if self.verify_spec():
                    self.uptodate = True
                    return
            except FileNotFoundError:
                pass
        self.uptodate = False
        self.data_dir = get_data_dir(self.path)
        self.get_last_activity()
        self.spec = {"data_dir": os.path.relpath(
            self.data_dir, start=self.path),
            "maxdim": "2d",
            "timestamp": self.timestamp}
        self.weakref = proxy(self)

        self._snapshots = None
        self._output_times = None
        self._snapshot_numbers = None
        self._output_times_dict = None
        self._fine_output_times = None

    def get(self, key=None, dim=None, planet=None, fluid=None):
        if dim is None:
            dim = self.spec["maxdim"]
        if fluid is None:
            fluid = "gas"

        if dim == "scalar":
            spec = self.spec["fluids"][fluid][dim][key]
            loader = loadscalar.ScalarLoader(
                spec["varname"], spec["datafile"], self.weakref)
            return loader()

    def scout(self):
        self.get_snapshots()
        self.get_units()
        self.get_domain_size()
        self.get_parameters()
        self.load_grid()
        self.load_times()
        self.get_fluids()
        self.get_planets()
        self.get_particles()
        self.get_fields()
        self.get_scalars()
        self.register_alias()

    def get_last_activity(self):
        """ Get the last time the simulation was updated.

        Access stats of the Quantities.dat file for the modified date.

        Raises
        ------
        FileNotFoundError
            If the proxy file does not exist.
        """
        file_path = os.path.join(self.data_dir, "monitor/Quantities.dat")
        self.timestamp = os.path.getmtime(file_path)

    def verify_spec(self):
        """ Make sure the spec applies to the current simulation.

        Check for the last modified time and compare to timestamp in spec.

        Returns
        -------
        bool
            True when spec is up to date, False otherwise.
        """
        try:
            if self.timestamp > self.spec["timestamp"]:
                return False
            else:
                return True
        except KeyError:
            return False

    @property
    def first_snapshot_number(self):
        if self._snapshots is None:
            self.get_snapshots()
        return self._snapshot_numbers[0]
    
    @property
    def last_snapshot_number(self):
        if self._snapshots is None:
            self.get_snapshots()
        return self._snapshot_numbers[-1]

    @property
    def snapshots(self):
        if self._snapshot_numbers is None:
            self.get_snapshots()
        return self._snapshot_numbers
    
    def get_snapshots(self):
        snapshot_time_file = self.datadir_path("snapshots/timeSnapshot.dat")
        self._snapshots = np.genfromtxt(snapshot_time_file, usecols=0, dtype=int)
        try:
            len(self._snapshots)
        except TypeError:
            self._snapshots = np.array([self._snapshots])
        self._snapshot_numbers = self._snapshots

    def get_parameters(self):
        if "paramters" in self.spec:
            self.parameters = self.spec["parameters"].copy()
            return
        try:
            param_file = self.datadir_path("snapshots/0/config.yml")
            self.parameters = loadparams.get_parameters(param_file)
            self.spec["parameters"] = self.parameters.copy()
        except FileNotFoundError:
            pass

    def get_domain_size(self):
        if all([key in self.spec for key in ["Nr", "Nphi"]]):
            self.Nr = self.spec["Nr"]
            self.Nphi = self.spec["Nphi"]
            return
        self.Nr, self.Nphi = np.genfromtxt(os.path.join(
            self.data_dir, "dimensions.dat"),
            usecols=(4, 5),
            dtype=int)
        self.spec["Nr"] = self.Nr
        self.spec["Nphi"] = self.Nphi

    def get_units(self):
        if "units" in self.spec:
            self.units = self.spec["units"].copy()
            return
        with open(self.datadir_path('units.yml'), 'r') as f:
            unitsdefs = yaml.safe_load(f)

        self.units = {}
        for key, val in unitsdefs.items():
            self.units[key] = val["unit"]

        self.spec["units"] = self.units.copy()

    def load_grid(self):
        try:
            if all([key in self.spec["grid"] for key in ["r_i", "phi_i"]]):
                self.r_i = quantity_from_spec(self.spec["grid"]["r_i"])
                self.phi_i = quantity_from_spec(self.spec["grid"]["phi_i"])
                return
        except KeyError:
            pass

        self.r_i = np.genfromtxt(self.data_dir +
                                 "/used_rad.dat") * u.Unit(self.units["length"])
        self.phi_i = np.linspace(0, 2 * np.pi, self.Nphi + 1) * u.rad
        self.spec["grid"] = {
            "r_i": (self.r_i.value, self.r_i.unit.to_string()),
            "phi_i": (self.phi_i.value, self.phi_i.unit.to_string())
        }

    def load_times(self):
        if all([key in self.spec for key in ["output_times", "fine_output_times"]]):
            self.output_times = quantity_from_spec(self.spec["output_times"])
            self.fine_output_times = quantity_from_spec(
                self.spec["fine_output_times"])
            return

        self._output_times = loadscalar.load_text_data_file(
            self.datadir_path("snapshots/timeSnapshot.dat"), "time")
        # self._snapshot_numbers = loadscalar.load_text_data_file(
        #     self.datadir_path("snapshots/timeSnapshot.dat"), "time step")
        self._output_times_dict = {
            n: t for n, t in zip(self._snapshot_numbers, self._output_times)}
        self._output_times_dict["reference"] = self._output_times[0]
        self._fine_output_times = loadscalar.load_text_data_file(
            self.datadir_path("monitor/Quantities.dat"), "time")
        self.spec["output_times"] = (
            self.output_times.value, self.output_times.unit.to_string())
        self.spec["fine_output_times"] = (
            self.fine_output_times.value, self.fine_output_times.unit.to_string())

    @property
    def output_times(self):
        if self._output_times is None:
            self.load_times()
        return self._output_times
    
    @property
    def snapshot_times(self):
        if self._output_times is None:
            self.load_times()
        return self._output_times
    
    def snapshot_time(self, N):
        if self._output_times_dict is None:
            self.load_times()
        return self._output_times_dict[N]

    @property
    def fine_output_times(self):
        if self._fine_output_times is None:
            self.load_times()
        return self._fine_output_times

    def register_alias(self):
        for planet in self.planets:
            planet.alias.register_dict(defs.alias_particle)
        self.fluids["gas"].alias.register_dict(defs.alias_fields)
        self.fluids["gas"].alias.register_dict(defs.alias_reduced)

    def get_particles(self):
        if "particles" in self.spec:
            pspec = self.spec["particles"]
            self.particles = loadparticles.ParticleLoader(
                pspec["name"],
                pspec["datafile_pattern"],
                quantity_from_spec(pspec["times"]),
                pspec["timesteps"], self.weakref)
            return
        p = re.compile(r"([\d]+)")
        timesteps = []
        datafile_pattern = "snapshots/{}/particles.dat"

        # exit if there is no particles file in the last snapshot
        last_particle_file = datafile_pattern.format(self._snapshot_numbers[-1])
        if not os.path.exists(os.path.join(self.data_dir, last_particle_file)):
            return

        for s in os.listdir(os.path.join(self.data_dir, "snapshots")):
            m = re.match(p, s)
            if m:
                key = m.groups()[0]
                if key == "autosave":
                    continue
                if os.path.exists(os.path.join(self.data_dir, datafile_pattern).format(key)):
                    timesteps.append(int(key))
        if len(timesteps) == 0:
            return
        timesteps.sort()
        imin = timesteps[0]
        times = u.Quantity([self.output_times[i-imin] for i in timesteps])
        self.particles = loadparticles.ParticleLoader(
            "dust", datafile_pattern, times, timesteps, self.weakref)
        self.spec["particles"] = {
            "name": "dust",
            "datafile_pattern": datafile_pattern,
            "times": (times.value, times.unit.to_string()),
            "timesteps": timesteps
        }

    def get_planets(self):
        if "planets" in self.spec:
            for val in self.spec["planets"].values():
                planet = particles.Planet(val["name"], val["pid"])
                for varname in val["variables"]:
                    datafile = val["variables"][varname]["datafile"]
                    loader = loadscalar.ScalarLoader(varname, datafile, self.weakref)
                    planet.register_variable(varname, loader)
                self.planets.append(planet)
            return
        self.spec["planets"] = {}
        planet_ids = []
        p = re.compile(r"(bigplanet|planet|nbody)(\d).dat")
        basename = ""
        for s in os.listdir(os.path.join(self.data_dir, "monitor")):
            m = re.match(p, s)
            if m:
                basename = m.groups()[0]
                planet_ids.append(m.groups()[1])
        planet_ids.sort()
        # create planets
        self.planets = []
        for pid in planet_ids:
            self.planets.append(particles.Planet(str(pid), pid))
            self.spec["planets"][str(pid)] = {
                "pid": pid, "name": str(pid), "variables": {}}
        # add variables to planets
        for pid, planet in zip(planet_ids, self.planets):
            planet_variables = loadscalar.load_text_data_variables(
                self.datadir_path("monitor/{}{}.dat".format(basename, pid)))
            for varname in planet_variables:
                datafile = "monitor/{}{}.dat".format(basename, pid)
                loader = loadscalar.ScalarLoader(varname, datafile, self.weakref)
                planet.register_variable(varname, loader)
                self.spec["planets"][str(pid)]["variables"][varname] = {
                    "datafile": datafile}

    def get_fluids(self):
        if "fluids" in self.spec:
            for fname in self.spec["fluids"]:
                self.fluids[fname] = fluid.Fluid(fname)
            return

        self.fluids["gas"] = fluid.Fluid("gas")
        self.spec["fluids"] = {"gas": {}}

    def get_fields(self):
        self.get_fields_2d()
        self.get_fields_1d()

    def get_fields_2d(self):
        if "fluids" in self.spec and all(("2d" in fspec for fspec in self.spec["fluids"].values())):
            for name, fspec in self.spec["fluids"].items():
                fl = self.fluids[name]
                for varname, vspec in fspec["2d"].items():
                    loader = load2d.FieldLoader2d(varname, vspec["info"], self.weakref)
                    fl.register_variable(varname, "2d", loader)
            return
        self.spec["fluids"]["gas"]["2d"] = {}
        gas = self.fluids["gas"]
        try:
            with open(self.datadir_path("info2D.yml"), "r") as f:
                infodict = yaml.safe_load(f)
        except FileNotFoundError:
            raise FileNotFoundError("Could not find 'info2D.yml' in simulatiuon output ('{}')".format(self.data_dir))
        
        if infodict is not None:
            for varname, info in infodict.items():
                loader = load2d.FieldLoader2d(varname, info, self.weakref)
                gas.register_variable(varname, "2d", loader)
                self.spec["fluids"]["gas"]["2d"][varname] = {"name": varname, "info": info}

    def get_fields_1d(self):
        if "fluids" in self.spec and all(("1d" in fspec for fspec in self.spec["fluids"].values())):
            for name, fspec in self.spec["fluids"].items():
                fl = self.fluids[name]
                for varname, vspec in fspec["1d"].items():
                    loader = load1d.FieldLoader1d(varname, vspec["info"], self.weakref)
                    fl.register_variable(varname, "1d", loader)
            return
        self.spec["fluids"]["gas"]["1d"] = {}
        gas = self.fluids["gas"]
        files = os.listdir(self.data_dir)
        for n in range(len(self.planets)):
            if "monitor/gas1D_torque_planet{}_0.dat".format(n) in files:
                varname = "torque planet {}".format(n)
                path_pattern = "monitor/gas1D_torque_planet{}_{}.dat".format(n, "{}")
                info = {
                    "pattern": path_pattern
                }
                loader = load1d.FieldLoader1dTorq(varname, info, self.weakref)
                gas.register_variable(varname, "1d", loader)
        if "gasMassFlow1D.info" in files:
            varname = "mass flow"
            info = {}
            loader = load1d.FieldLoader1dMassFlow(varname, info, self.weakref)
            gas.register_variable(varname, "1d", loader)
        for fname in files:
            m = re.search(r"(.*)1D\.info", fname)
            if m:
                stem = m.groups()[0]
                infofile = fname
                varname = stem.lower()
                if varname == "massflow":
                    continue
                info = {"infofile": infofile}
                loader = load1d.FieldLoader1d(
                    varname, info, self.weakref)
                gas.register_variable(varname, "1d", loader)
                self.spec["fluids"]["gas"]["1d"][varname] = {
                    "varname": varname, "info": info}

    def get_scalars(self):
        if "fluids" in self.spec and all(("scalar" in fspec for fspec in self.spec["fluids"].values())):
            for name, fspec in self.spec["fluids"].items():
                fl = self.fluids[name]
                if "scalar" in fspec:
                    sspec = fspec["scalar"]
                    for varname, vspec in sspec.items():
                        loader = loadscalar.ScalarLoader(
                            varname, vspec["datafile"], self.weakref)
                        fl.register_variable(varname, "scalar", loader)
            return

        self.spec["fluids"]["gas"]["scalar"] = {}
        gas = self.fluids["gas"]
        for datafile in ["monitor/Quantities.dat", "monitor/timestepLogging.dat"]:
            variables = loadscalar.load_text_data_variables(
                self.datadir_path(datafile))
            for varname, _ in variables.items():
                loader = loadscalar.ScalarLoader(varname, datafile, self.weakref)
                gas.register_variable(varname, "scalar", loader)
                self.spec["fluids"]["gas"]["scalar"][varname] = {
                    "varname": varname, "datafile": datafile}