from collections import OrderedDict
import importlib
import sys

class UnknownCodeError(Exception):
    pass

class MultipleCodeError(Exception):
    pass

# import all loaders
import loader

def identify_code(path):
    code_list = []
    for key, mod in loader.available.values():
        if mod.identify_code(path):
            code_list.append(code)
    if len(code_list) == 0:
        raise UnknownCodeError("No known code matches data in '{}'".format(path))
    elif len(code_list) > 1:
        raise MultipleCodeError("Multiple codes identified the data in '{}' which where '{}'".format(path, code_list))
    return code_list[0]

class Simdata:
    def __init__(self, path, loader=None, **kwargs):
        self.path = path
        if loader:
            if "Loader" in dir(loader):
                # assume its a module
                self.code = loader.code_info
                self.loader = loader.Loader(path, **kwargs)
            else:
                # assume its a loader object
                self.loader = loader
                self.code = type(loader).code_info
        else:
            self.code = self.identify_code()
            self.loader = loader.available[code].Loader(path, **kwargs)
        self.fields = {}
        self.nbodysystems = {}
        self.parameters = {}
        self.meta = {}

    
