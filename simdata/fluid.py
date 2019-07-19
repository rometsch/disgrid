from .alias import Alias
supported_variable_types = ["field", "vector"]

class Fluid:
    def __init__(self, name, alias = Alias()):
        self.name = name
        self.variable_loaders = {"field": {}, "vector": {}}
        self.alias = alias

    def register_variable(self, name, variable_type, loader_function):
        if not variable_type in supported_variable_types:
            raise ValueError("Variable type '{}' for '{}' not supported".format(variable_type, name))

        if not callable(loader_function):
            raise TypeError("Loader function '{}' for '{}' not callable".format(loader_function, name))

        self.variable_loaders[variable_type][name] = loader_function

    def register_alias(self, new_alias):
        for key, item in new_alias.items():
            self.alias.register(key, item)

    def get(self, variable_type, name, num_output=None, *args, **kwargs):
        name = self.alias(name)
        if variable_type == "field":
            if num_output is None:
                raise TypeError("get() missing 1 required optional argument: 'num_output' for variable_type=field")
        if variable_type == "field":
            return self.variable_loaders[variable_type][name](num_output, *args, **kwargs)
        elif variable_type == "vector":
            return self.variable_loaders[variable_type][name](*args, **kwargs)
        else:
            raise TypeError("Unknown variable_type '{}'".format(variable_type))
