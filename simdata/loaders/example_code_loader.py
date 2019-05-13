code_info = ( "exmaple_code", "0.1", "test")

import os
from . import interface

def identify(path):
    try:
        with open(os.path.join(path, "text.txt")) as f:
            for line in f:
                if line.strip() == "very specific content!":
                    return True
            return False
    except FileNotFoundError:
        return False

class Loader(interface.LoaderInterface):

    def __init__(self, path):
        self.grid = None

    def load_field_rho_data(self, n):
        return np.fromfile(os.path.join(self.dataDir, "gasdens{}.dat".format(n)))

    def load_field_rho_time(self, n):
        return self.times[n]

    def load_field_rho_grid(self, n):
        return self.grid

    def load_grid(self):
        pass

class FieldRho(interface.FieldQuantity):

    def load_data(self, n):
        pass
