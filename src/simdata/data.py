from collections import OrderedDict
import importlib
import sys
import os
import copy

# import all loaders
from . import loaders


class Data:
    """ Create a data interface for the data in 'path' """
    def __init__(self, path=os.getcwd(), loader=None, **kwargs):
        self.code, self.loader = loaders.get_loader(path, loader, **kwargs)
        self.loader.scout()
        self.fluids = self.loader.fluids
        self.particlegroups = self.loader.particlegroups
        self.planets = self.loader.planets
        self.parameters = self.loader.parameters

    def get_fluid(self, name):
        return self.fluids[name]
