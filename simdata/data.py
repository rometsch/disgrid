from collections import OrderedDict
import importlib
import sys
import os
import copy

# import all loaders
from . import loaders

class Data:
    def __init__(self, path=os.getcwd(), data_dir=None, loader=None, search_args=None, **kwargs):
        self.path = try_simscripts_lookup(path, search_args)
        self.code, self.loader = loaders.get_loader(self.path, loader, **kwargs)
        if data_dir is not None:
            self.loader.set_data_dir(data_dir)
        self.loader.scout()
        self.fluids = self.loader.fluids
        self.particlegroups = self.loader.particlegroups
        self.planets = self.loader.planets
        self.parameters = self.loader.parameters
        self.meta = self.loader.meta

    def get_fluid(self, name):
        return self.fluids[name]

def try_simscripts_lookup(pattern, search_args=None):
    try:
        import simscripts.cache
        c = simscripts.cache.Cache()
        try:
            if search_args is not None:
                rv = c.search(pattern, **search_args)["path"]
            else:
                rv = c.search(pattern)["path"]
        except simscripts.cache.NoSimulationFoundError:
            rv = pattern
        except simscripts.cache.ResultNotUniqueError:
            raise
    except ImportError:
        rv = pattern
    return rv
