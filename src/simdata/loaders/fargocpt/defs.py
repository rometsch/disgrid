""" Definitions for the simdata fargocpt loader.

This file contains dicts defining the output variables of fargocpt, their naming patterns and units.
"""

vars2d = {
    "dens": {
        "pattern": "gasdens{}.dat",
        "unit": "g cm-2"
    },
    "energy": {
        "pattern": "gasenergy{}.dat",
        "unit": "erg"
    },
    "vrad": {
        "pattern": "gasvrad{}.dat",
        "unit": "cm s-1",
        "interfaces": ["r"]
    },
    "vtheta": {
        "pattern": "gasvtheta{}.dat",
        "unit": "cm s-1",
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
