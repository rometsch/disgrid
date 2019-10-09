from .alias import Alias
supported_geometries= ["3d", "2d", "1d", "scalar"]

class Fluid:
    def __init__(self, name, alias = Alias()):
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
        if geometry == "scalar":
            return self.variable_loaders[geometry][name](*args, **kwargs)
        elif geometry in supported_geometries:
            if num_output is None:
                raise TypeError("get() missing 1 required optional argument: 'num_output' for geometry={}".format(", ".join(supported_geometries[:-1])))
            return self.variable_loaders[geometry][name](num_output, *args, **kwargs)
        else:
            raise TypeError("Unknown geometry '{}'".format(geometry))
