from collections import OrderedDict
import importlib
import sys

class UnknownCodeError(Exception):
    pass

class MultipleCodeError(Exception):
    pass

# import all loaders
from . import loaders

def identify_code(path):
    code_list = []
    for key, mod in loaders.available.items():
        if mod.identify(path):
            code_list.append(key)
    if len(code_list) == 0:
        raise UnknownCodeError("No known code matches data in '{}'".format(path))
    elif len(code_list) > 1:
        raise MultipleCodeError("Multiple codes identified the data in '{}' which where '{}'".format(path, code_list))
    return code_list[0]

def get_loader(path, loader, **kwargs):
    if loader:
        if "Loader" in dir(loader):
            # assume its a module
            code = loader.code_info
            loader = loader.Loader(path, **kwargs)
        else:
            # assume its a loader object
            loader = loader
            code = type(loader).code_info
    else:
        code = identify_code(path)
        loader = loaders.available[code].Loader(path, **kwargs)
    return code, loader



class Data:
    def __init__(self, path, loader=None, **kwargs):
        self.path = path
        self.code, self.loader = get_loader(path, loader, **kwargs)
        self.fluids = {}
        self.particlegroups = {}
        self.planets = {}
        self.parameters = {}
        self.meta = {}
        self.grids = {}

    def get_fluid(name):
        return self.fluids[name]


    
