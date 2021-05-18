# a loader object orchestrates the loading of data from files
# it
# it provides a framework for caching data and provides function which return field objects
import os
import shutil

import numpy as np

from .. import field
from ..filecache import FileCache


class Interface:
    def __init__(self, path, owner=None, file_caching=True, update=True, **kwargs):
        self.path = path
        self.fluids = {}
        self.scalar = {}
        self.particles = {}
        self.particlegroups = {}
        self.planets = []
        self.parameters = {}
        
        # caching
        self.owner = owner
        self.file_caching = file_caching
        self.update = update
        self.init_cache(owner, file_caching, update)

    def scout(self):
        # find all variables
        # adjust self.fluids
        # adjust self.field_loaders
        pass

    def get(self, *args, **kwargs):
        pass

    def get_output_time(self, n):
        pass

    def cached(self, filename, changing=False):
        if not self.file_caching:
            return os.path.join(self.data_dir, filename)

        if os.path.isabs(filename):
            if os.path.commonpath([os.path.abspath(self.data_dir), os.path.abspath(filename)]) == self.data_dir:
                filename = os.path.relpath(filename, self.data_dir)

        # don't reuse an old changing file
        if changing and not hasattr(self, "uptodate") and not filename in self.cached_files:
            invalidate_entry = True
        # be aware of the uptodate flag
        elif changing and not self.uptodate and not filename in self.cached_files:
            invalidate_entry = True
        else:
            invalidate_entry = False
        
        if not self.update:
            invalidate_entry = False

        filepath = self.filecache.cached_file(filename, self.data_dir, invalidate=invalidate_entry)
        
        self.cached_files.add(filename)
        
        return filepath

    def init_cache(self, owner, file_caching, update):
        if not self.file_caching:
            self.cached = self.datadir_path
        else:
            try:
                simid = self.owner.sim["uuid"]
                if os.path.exists(self.owner.sim["path"]):
                    raise AttributeError()  # its a local path so step out of try
                cachedir_base = self.owner.config["cachedir"]
                self.cachedir = os.path.join(cachedir_base, simid)
                self.filecache = FileCache(self.cachedir, simid)
            except (KeyError, AttributeError):
                self.cached = self.datadir_path
        
        self.owner = owner
        self.file_caching = file_caching
        self.cached_files = set()

    def datadir_path(self, filename, **kwargs):
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
