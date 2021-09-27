import os

# import all loaders
from . import loaders


class Data:
    """ Create a data interface for the data in 'path' """

    def __init__(self, path=os.getcwd(), loader=None, init_hooks=None, **kwargs):
        self.path = path
        self.loader = loader
        self.kwargs = kwargs
        self.init_hooks = init_hooks

    def init(self):
        self.code, self.loader = loaders.get_loader(
            self.path, self.loader, **self.kwargs)
        self.loader.scout()
        self.fluids = self.loader.fluids
        self.particles = self.loader.particles
        self.particlegroups = {"particles": self.particles}
        self.planets = self.loader.planets
        self.parameters = self.loader.parameters
        self.register_postprocessor()
        if self.init_hooks is not None:
            for func in self.init_hooks:
                func()

    def __getattr__(self, attr):
        try:
            return self.__dict__[attr]
        except KeyError:
            if attr in ["path", "code", "loader", "fluids",
                        "particles", "particlegroups",
                        "planets", "paramters"]:
                self.init()
                return self.__dict__[attr]
            else:
                raise AttributeError(f'{attr}')

    def get_fluid(self, name):
        return self.fluids[name]

    def register_postprocessor(self):
        """ Add new loaders to calculate derived quantities. """
        try:
            from simprocess import register
            register(self)
        except ImportError:
            pass

    def get(self, var=None, dim=None, N=None, planet=None, t=None, fluid="gas", **kwargs):
        """ General getter function. """
        if planet is not None:
            return self.planets[planet].get(var, **kwargs)
        else:
            return self.fluids[fluid].get(dim, var, num_output=N, t=t, **kwargs)