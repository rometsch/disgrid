""" Functions to load parameters from fargo3d output files using the summmary files.
"""
import yaml


def get_parameters(param_file):
    """ Read simulation parameters from a fargocpt parameter file.

    Fargocpt uses an ini style parameter file.

    Parameters
    ----------
    param_file: str
        Path to the parameter file.

    Returns
    -------
    dict
        Dictionary contining parameters as key value pairs.
    """
    with open(param_file) as f:
       rv = yaml.safe_load(f)
    return rv
