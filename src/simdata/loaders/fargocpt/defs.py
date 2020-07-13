""" Definitions for the simdata fargocpt loader.

This file contains dicts defining the output variables of fargocpt, their naming patterns and units.
"""

import astropy.units as u

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
    "vtheta": {
        "pattern": "gasvtheta{}.dat",
        "unit": u.cm / u.s,
        "interfaces": ["phi"]
    }
}

alias_fields = {
    "mass density": "dens",
    "velocity radial": "vrad",
    "velocity azimuthal": "vtheta",
    "energy density": "energy"
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
