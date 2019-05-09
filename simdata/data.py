from collections import OrderedDict
import importlib
import sys

class UnknownCodeError(Exception):
    pass

class MultipleCodeError(Exception):
    pass

# import all loaders
from . import loaders

class Data:
    def __init__(self, path, loader=None, **kwargs):
        self.path = path
        self.code, self.loader = loaders.get_loader(path, loader, **kwargs)
        self.fluids = {}
        self.particlegroups = {}
        self.planets = {}
        self.parameters = {}
        self.meta = {}
        self.grids = {}

    def get_fluid(name):
        return self.fluids[name]


    
