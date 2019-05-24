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
        self.loader.scout()
        self.fluids = self.loader.fluids
        self.vectors = self.loader.vectors
        self.particles = self.loader.particles
        self.planets = self.loader.planets
        self.parameters = self.loader.parameters
        self.meta = self.loader.meta

    def get_fluid(name):
        return self.fluids[name]


    
