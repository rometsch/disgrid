""" Functions to load particle data from fargocpt output files.
"""
import astropy.units as u
import numpy as np


def load_particle_data(datafile, units):
    """ Load a variable from a fargocpt datafile for all particles.

    Parameters
    ----------
    datafile : str
        Path to the datafile.

    varname : str
        Name of the variable.

    Returns
    -------
    dict
        Dictionary containing id, x, y and size (the particle radius).
    """
    res = np.fromfile(
        datafile, dtype=[('Id', np.dtype(int)), ('Values', np.dtype(float), 11)])
    ids = res["Id"]
    vals = res["Values"]

    particles = {
        "id": ids,
        "r": vals[:, 0]*u.Unit(units["length"]),
        "phi": vals[:, 1]*u.rad,
        "r dot": vals[:, 2]*u.Unit(units["length"])/u.Unit(units["time"]),
        "phi dot": vals[:, 3]*u.rad/u.Unit(units["time"]),
        "r ddot": vals[:, 4]*u.Unit(units["length"])/u.Unit(units["time"])**2,
        "phi ddot": vals[:, 5]*u.rad/u.Unit(units["time"])**2,
        "mass": vals[:, 6]*u.Unit(units["mass"]),
        "size": vals[:, 7]*u.Unit(units["length"]),
        "timestep": vals[:, 8],
        "facold": vals[:, 9],
        "stokes": vals[:, 10]
    }
    return particles


class ParticleLoader:
    def __init__(self, name, datafile_pattern, time, timesteps, loader, *args, **kwargs):
        self.loader = loader
        self.datafile_pattern = datafile_pattern
        self.time = time
        self.timesteps = timesteps
        self.time_dict = {N: t for N, t in zip(timesteps, time)}
        self.name = name

    def get(self, N):
        time = self.load_time(N)
        data = self.load_data(N)
        return data, time

    def load_data(self, N):
        filepath = self.loader.datadir_path(self.datafile_pattern.format(N))
        rv = load_particle_data(
            filepath, self.loader.units)
        return rv

    def load_time(self, N):
        rv = self.time_dict[N]
        return rv
