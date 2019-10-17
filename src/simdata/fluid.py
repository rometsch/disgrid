from .alias import Alias
supported_geometries= ["3d", "2d", "1d", "scalar"]

class Fluid:
    def __init__(self, name, alias=None):
        if alias is None:
            alias = Alias()
        self.name = name
        self.variable_loaders = {g:{} for g in supported_geometries}
        self.alias = alias

    def register_alias(self, new_alias):
        for key, item in new_alias.items():
            self.alias.register(key, item)

    def register_variable(self, name, geometry, loader_function):
        if not callable(loader_function):
            raise TypeError("Loader function '{}' for '{}' not callable".format(loader_function, name))
        if not geometry in supported_geometries:
            raise ValueError("Geometry '{}' for '{}' not supported".format(geometry, name))
        self.variable_loaders[geometry][name] = loader_function

    def get(self, geometry, name, num_output=None, *args, **kwargs):
        name = self.alias(name)
        loader = self._get_loader(geometry, name)
        if geometry == "scalar":
            return loader(*args, **kwargs)
        elif geometry in supported_geometries:
            if num_output is None:
                raise TypeError("get() missing 1 required optional argument: 'num_output' for geometry={}".format(", ".join(supported_geometries[:-1])))
            return loader(num_output, *args, **kwargs)

    def _get_loader(self, geometry, name):
        if geometry in supported_geometries:
            return self.variable_loaders[geometry][name]
        else:
            raise TypeError("Unknown geometry '{}'".format(geometry))
        
    def get_time(self, geometry, name, num_output=None):
        """Return an array containing the time for each output for the given quantity"""
        name = self.alias(name)
        loader = self._get_loader(geometry, name)
        return loader.load_time(num_output)
        
    def get_grid(self, geometry, name, num_output=0):
        """Get a grid appropriate for the given quantity"""
        name = self.alias(name)
        loader = self._get_loader(geometry, name)
        return loader.load_grid(num_output)

