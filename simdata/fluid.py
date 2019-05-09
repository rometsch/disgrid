supported_variable_types = ["field", "reduced"]


class Fluid:
    def __init__(self, name, variable_loaders={"field": {}, "reduced": {}}):
        self.name = name
        if not all([key in supported_variable_types for key in variable_loaders]):
            raise ValueError("One of types '{}' not supported".format([key for key in variable_loaders]))
        self.variable_loaders = variable_loaders

    def register_variable(self, name, variable_type, loader_function):
        if not variable_type in supported_variable_types:
            raise ValueError("Variable type '{}' for '{}' not supported".format(variable_type, name))
        
        if not callable(loader_function):
            raise TypeError("Loader function '{}' for '{}' not callable".format(loader_function, name))

        self.variable_loaders[variable_type][name] = loader_function

    def get(self, variable_type, name, *args, **kwargs):
        return self.variable_loaders[varialbe_type][name](*arg, **kwargs)

    def get_time(self, variable_type, name):
        return self.variable_loader[variable_type][name](*args, **kwargs, return_times=True)



