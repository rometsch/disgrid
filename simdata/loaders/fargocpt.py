code_info = ( "fargocpt", "0.1", "testloader")

import os
from . import interface

def identify(path):
    return "fargo" in os.listdir(path)

vars2d = {
    "gasdens" : "gasdens{}.dat",
    "gasenergy" : "gasenergy{}.dat"
    }

vars1d = {
    "mass" : ""
    }

def var_in_files(var, files):
    return True

def load_scalar(file, var):
    return [1,1]

class Loader(interface.Interface):
    def scout(self):
        self.fluids["gas"] = fluid.Fluid("gas")
        gas = self.fluids["gas"]
        files = os.listdir(self.path)
        for v, p in vars2d.items():
            if var_in_files(v, files):
                gas.register_variable(v, "field", FieldLoader(v, d, self))



class FieldLoader:

    def __init__(self, name, pattern, loader, *args, **kwargs):
        self.loader = loader
        self.pattern = pattern
        self.name = name

    def __call__(self, n):
        t = self.loader.get_time(n)
        f = field.Field(self.load_grid(), self.load_data(), t, self.name)
        return f

    def load_data(self, n):
        dataDir = self.loader.dataDir
        rv = np.fromfile(dataDir + "/" + pattern.format(self.n)).reshape(self.loader.Nr, self.loader.Nphi)
        return rv

    def load_grid(self, n):
        r_i = np.genfromtxt(dataDir + "/used_rad.dat")
        Nphi = 101
        phi_i = np.linspace(-np.pi, np.pi, Nphi+1)
        g = grid.Grid(r_i = r_i, phi_i = phi_i)
        return g
