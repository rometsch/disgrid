""" Definitions for the example loader.

This file contains dicts defining the output variables in the example dataset, their naming patterns and units.
"""

vars_2d = {
    'density': {
        'file': 'density{}.dat',
        'unitpowers': {
            "mass": 1,
            "length": -2
        }
    },
    'vx': {
        'file': 'vx{}.dat',
        'unitpowers': {
            "time": -1,
            "length": 1
        },
        'interfaces' : ['x']
    },
    'vy': {
        'file': 'vy{}.dat',
        'unitpowers': {
            "time": -1,
            "length": 1
        },
        'interfaces' : ['y']
    },
}

vars_scalar = {
    'mass': {
        'file': 'quantities.dat',
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {
            "mass": 1
        }
    },
    'momentum': {
        'file': 'quantities.dat',
        'datacol': 2,
        'timecol': 0,
        'unitpowers': {
            "length": 1,
            "time": -1
        }
    }
}

planet_vars_scalar = {
    'x': {
        'file': 'planet.dat',
        'datacol': 1,
        'timecol': 0,
        'unitpowers': {
            'length': 1
        }
    },
    'y': {
        'file': 'planet.dat',
        'datacol': 2,
        'timecol': 0,
        'unitpowers': {
            'length': 1
        }
    }
}
