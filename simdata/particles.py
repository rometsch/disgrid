from collections import OrderedDict
from . import vector
from .aliases import Aliases

required_data = ["time", "position", "velocity"]

def check_data(data):
    # need at least time, position and velocities
    for varname in required_data:
        if not varname in data:
            raise KeyError("nbody data must contain "+varname)

class NbodySystem:
    # class to hold data for an nbody system for all timesteps
    def __init__(self, name):
        self.name = name
        self.ids = []
        self.variable_loaders = {}
        self.aliases = Aliases()

    def register_particles(self, ids):
        for i in ids:
            self.ids.append(i)

    def register_alias(self, new_aliases):
        for key, item in new_aliases.items():
            self.aliases.register(key, item)

    def register_variable(self, name, loader_function):
        if not callable(loader_function):
            raise TypeError("Loader function '{}' for '{}' not callable".format(loader_function, name))

        self.variable_loaders[name] = loader_function

    def get(self, name, num_particles=None, *args, num_output=None, **kwargs):
        name = self.aliases(name)
        if num_particles is None:
            num_particles = slice(len(self.ids))
        return self.variable_loaders[name](num_output, self.ids[num_particles], *args, **kwargs)
        
class ParticleVector(vector.Vector):
    pass
    
    # def get(self, ids=None, *args, **kwargs):
    #     data = super().get(*args, **kwargs)
    #     time = None
    #     if isinstance(data, tuple):
    #         time = data[0]
    #         data = data[1]
    #     # select particles
    #     if len(data.shape) == 2:
    #         rv = data[:,ids]
    #     else:
    #         rv = data[:,:,ids]
    #     if time is not None:
    #         return (time, rv)
    #     else:
    #         return rv
