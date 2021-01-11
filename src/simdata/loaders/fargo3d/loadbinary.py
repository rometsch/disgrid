""" Load binary output data for FARGO3D.

FARGO3D supports different output formats of binary data.

The legacy output format uses separate files for each field.
The mpiio format saves all fields into a single file.

"""
import os

import numpy as np

from .mpiio import Fields


def load_data(filename, varname, n):
    """ Load binary data of a variable at step n.

    Output type, legacy or mpiio, is automatically determined using the file extension.

    Parameteters
    ------------
    filename : str
        Path to the data file.
    varname : str
        Name of the variable (needed for mpiio format).
    n : int
        Output number.

    Returns
    -------
    np.array
        1D flat output array.
    """
    file_ext = filename.split(".")[-1]
    if file_ext == "mpiio":
        print(f"loading {varname} from mpiio datafile {filename} at output {n}")
        directory = os.path.dirname(filename)
        fluidname = os.path.basename(filename).split("_")[0]
        data = Fields(directory, fluidname, n).get_field(varname)
    else:
        data = np.fromfile(filename)
    return data
