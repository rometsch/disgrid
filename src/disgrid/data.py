import os
from disgrid.loaders import get_loader
from copy import deepcopy

class Data:
    """ Create a data interface for the data in 'path' """

    def __init__(self, path=os.getcwd(), loader_hint=None, init_hooks=None, **kwargs):
        if not os.path.exists(path):
            raise FileNotFoundError(f"The path {path} does not exist!")
        self.path = path
        self._loader_hint = loader_hint
        self._loader = None
        self.kwargs = kwargs
        self.init_hooks = init_hooks

    def init(self):
        # import all loaders
        _, self._loader = get_loader(
            self.path, self._loader_hint, **self.kwargs)
        self.loader.scout()
        self.register_postprocessor()
        if self.init_hooks is not None:
            for func in self.init_hooks:
                func()

    @property
    def loader(self):
        if self._loader is None:
            self.init()
        return self._loader

    @property
    def code(self):
        return deepcopy(self.loader.code)
    
    @property
    def fluids(self):
        return self.loader.fluids
    
    @property
    def particles(self):
        return self.loader.particles
    
    @property
    def particlegroups(self):
        return {"particles": self.loader.particles}
    
    @property
    def planets(self):
        return self.loader.planets
    
    @property
    def parameters(self):
        return self.loader.parameters

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
            return self.planets[int(planet)].get(var, **kwargs)
        else:
            # select highest dimension if none is given
            if dim is None:
                dim = [key for key, val in self.fluids[fluid].variable_loaders.items() if len(val) > 0][0]
            return self.fluids[fluid].get(dim, var, num_output=N, t=t, **kwargs)
        
    def avail(self):
        """ Tell what data is available. """

        fluids = {}
        for fluid in self.fluids:
            fluids[fluid] = {}
            for k,v in self.fluids[fluid].variable_loaders.items():
                names = [k for k in v.keys()]
                if len(names) > 0:
                    fluids[fluid][k] = names

        planets = {}
        for n,planet in enumerate(self.planets):
            names = [k for k in planet.variable_loaders.keys()]
            planets[f"{n}"] = names
        
        Nsnapshots = len(self.loader.output_times)
        rv = {
            "fluids" : fluids,
            "planets" : planets,
            "Nsnapshots" : Nsnapshots,
            "code" : deepcopy(self.loader.code_info)
        }
        try:
            rv["Nfirst"] = self.loader.first_snapshot
            rv["Nlast"] = rv["Nfirst"] + Nsnapshots -1
        except AttributeError:
            rv["Nfirst"] = 0
            rv["Nlast"] = Nsnapshots -1
        
        return rv