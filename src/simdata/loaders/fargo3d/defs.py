""" Definitions for the simdata fargo3d loader.

This file contains dicts defining the output variables of fargo3d, their naming patterns and units.
"""
vars_maxdim = {
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
    "velocity radial": {
        "pattern": "{}vy{}.dat",
        "unitpowers": {
            "length": 1,
            "time": -1
        },
        "interfaces": ["r"],
    },
    "velocity azimuthal": {
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
        'pattern': '{}vy{}.dat',
        'directions': ["r"],
        'unitpowers': {
            "mass": 0,
            "length": 1,
            "time": -1
        }
    },
    'velocity azimuthal': {
        'pattern': '{}vx{}.dat',
        'directions': ["r"],
        'unitpowers': {
            "mass": 0,
            "length": 1,
            "time": -1
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
    'omega frame': {
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
    # orbital elements
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
