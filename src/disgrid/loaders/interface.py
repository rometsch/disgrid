# a loader object orchestrates the loading of data from files
# it
# it provides a framework for caching data and provides function which return field objects
import os
import shutil

import numpy as np

from .. import field

class Interface:
    def __init__(self, path, update=True, **kwargs):
        self.path = path
        self.fluids = {}
        self.scalar = {}
        self.particles = {}
        self.particlegroups = {}
        self.planets = []
        self.parameters = {}
        
        self.update = update

    def scout(self):
        # find all variables
        # adjust self.fluids
        # adjust self.field_loaders
        pass

    def get(self, *args, **kwargs):
        pass

    def get_output_time(self, n):
        pass

    def datadir_path(self, filename):
        """ Return path of the file inside the data directory.
        
        Used for the case with caching disabled.
        
        Parameteters
        ------------
        str
            Name of the datafile.        

        Returns
        -------
        str
            Filepath inside the datadir.
        """
        filepath = os.path.join(self.data_dir, filename)
        return filepath


class FieldLoader:
    def __init__(self, name, info, loader, *args, **kwargs):
        self.loader = loader
        self.info = info
        self.name = name

    def __call__(self, n, *args, **kwargs):
        f = field.Field(self.load_grid(n), self.load_data(n),
                        self.load_time(n, *args, **kwargs), self.name)
        return f

    def load_time(self, n):
        raise NotImplementedError(
            "This is a virtual method. Please use a FieldLoader{}d for the specific geometry"
        )

    # def load_times(self,):
    #     """Returns an array containing the time for each output"""
    #     return self.load_time(slice(0,-1,1))

    def load_data(self, n):
        raise NotImplementedError(
            "This is a virtual method. Please use a FieldLoader{}d for the specific geometry"
        )

    def load_grid(self, n):
        raise NotImplementedError(
            "This is a virtual method. Please use a FieldLoader{}d for the specific geometry"
        )
