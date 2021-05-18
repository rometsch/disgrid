# a loader object orchestrates the loading of data from files
# it
# it provides a framework for caching data and provides function which return field objects
import os
import shutil

import numpy as np

from .. import field


class Interface:
    def __init__(self, path, owner=None, file_caching=False, **kwargs):
        self.path = path
        self.fluids = {}
        self.scalar = {}
        self.particles = {}
        self.particlegroups = {}
        self.planets = []
        self.parameters = {}
        self.owner = owner
        self.file_caching = file_caching
        self.cached_files = []

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
            use_cache = False
        # be aware of the uptodate flag
        elif changing and not self.uptodate and not filename in self.cached_files:
            use_cache = False
        else:
            use_cache = True

        data_dir = self.data_dir

        if use_cache:
            try:
                simid = self.owner.sim["uuid"]
                if os.path.exists(self.owner.sim["path"]):
                    raise AttributeError()  # its a local path so step out of try
                cachedir_base = self.owner.config["cachedir"]
                cachedir = os.path.join(cachedir_base, simid)
                os.makedirs(cachedir, exist_ok=True)
                filepath_in_src = os.path.join(self.data_dir, filename)
                old_filename = filename
                if ".." in filename:
                    filename = filename.replace("..", "__subdir__")
                filepath_in_cache = os.path.join(cachedir, filename)
                if not os.path.exists(filepath_in_cache):
                    os.makedirs(os.path.dirname(
                        filepath_in_cache), exist_ok=True)
                    shutil.copy2(filepath_in_src, filepath_in_cache)
                data_dir = cachedir
                self.cached_files.append(old_filename)
            except (KeyError, AttributeError):
                pass
        
        filepath = os.path.join(data_dir, filename)
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
