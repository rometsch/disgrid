""" Functions to load parameters from fargo3d output files using the summmary files.
"""


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
    parameters = {}
    with open(param_file) as f:
        for line in f:
            line = line.strip()
            if line == "" or line[0] in ["#", "="]:
                continue
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
