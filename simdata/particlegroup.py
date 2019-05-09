from collections import OrderedDict

required_data = ["time", "x", "y"]

def check_data(self):
    # need at least time, position and velocities
    for varname in required_data:
        if not varname in self.data:
            raise KeyError("nbody data must contain "+varname)

class NbodySystem:
    # class to hold data for an nbody system at one timestep
    def __init__(self, dim):
        self.dim = dim
        self.data = OrderedDict() # array of dim x Nparticles

class NbodySystemFull:
    # class to hold data for an nbody system for all timesteps
    def __init__(self, dim):
        self.dim = dim
        self.data = OrderedDict() # array of dim x Nparticles x Ntimesteps

    def load(self, varname, n, tmax=None):
        if tmax is not None:
            tmin = n
        elif n.size == 2:
            tmin, tmax = n

        if tmax:
            self.load_range(varname, tlim)
        else:
            self.load_single(varname, n)

    def load_range(self, varname, tlim):
        raise NotImplementedError("load_range function must be implemented in code specific loader")

    def load_range(self, varname, n):
        raise NotImplementedError("load_sing function must be implemented in code specific loader")
