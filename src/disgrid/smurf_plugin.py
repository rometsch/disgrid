# Integration of disgrid into smurf.
# disgrid is a python package to handle simulation data.
# This module provides the functionality to use sim ids
# to transparently get a data object without manually
# looking up simulation dir paths.

import os
import pickle
import time
import numpy as np

from pathlib import Path

from . import data
from .config import Config
from smurf.cache import LocalSimCache, get_cache_by_id, CacheMiss
from smurf.info import Info
from smurf.search import search

class SmurfData(data.Data):
    # data loader with smurf support to locate remote simulations
    # and mount them via sshfs
    def __init__(self, simid, search_remote=True, search_args={}, **kwargs):
        simid = insert_local_sim_to_cache(simid)
        self.sim = search(simid,
                          remote=search_remote,
                          unique=True,
                          **search_args)[0]
        path = self.sim["path"]
        if "disgrid_code" in self.sim:
            kwargs["loader"] = self.sim["disgrid_code"]
        elif "simcode" in self.sim:
            kwargs["loader"] = self.sim["simcode"]
        self.simid = simid
        self.config = Config()
        if "spec" in kwargs:
            spec = kwargs["spec"]
        else:
            spec = self.search_spec()
        if spec is not None:
            kwargs["spec"] = spec
        super().__init__(path, init_hooks=[
            self.register_code, self.save_spec], owner=self, file_caching=True, **kwargs)

    def register_code(self):
        self.sim["disgrid_code"] = self.code
        c = get_cache_by_id(self.sim["uuid"])
        if c is not None:
            c.insert(self.sim["uuid"], self.sim)

    def search_spec(self):
        spec = None
        simid = self.sim["uuid"]
        if "specdir" in self.config.data:
            specdir = self.config.data["specdir"]
        elif "cachedir" in self.config.data:
            specdir = self.config.data["cachedir"]
        else:
            return None
        specfile = os.path.join(specdir, f"{simid}.spec.pickle")
        if os.path.exists(specfile):
            try:
                spec = load_pickle(specfile)
                if "timestamp" in spec:
                    self._old_spec_timestamp = spec["timestamp"]
            except EOFError:
                pass
        return spec

    def save_spec(self):
        if hasattr(self, "_old_spec_timestamp") and hasattr(self.loader, "spec"):
            if "timestamp" in self.loader.spec:
                if self.loader.spec["timestamp"] <= self._old_spec_timestamp:
                    return
        simid = self.sim["uuid"]
        if "specdir" in self.config.data:
            specdir = self.config.data["specdir"]
        elif "cachedir" in self.config.data:
            specdir = self.config.data["cachedir"]
        else:
            return
        if hasattr(self.loader, "spec"):
            specfile = os.path.join(specdir, f"{simid}.spec.pickle")
            save_pickle(specfile, self.loader.spec)


def insert_local_sim_to_cache(path):
    """ Try to insert a local directory to the cache if not present. """
    abspath = os.path.abspath(path)
    if os.path.isdir(os.path.abspath(path)):
        try:
            info = Info(abspath)
            try:
                lc = LocalSimCache()
                lc.request(info.uuid)
                rv = info.uuid
            except CacheMiss:
                rv = lc.add_sim_to_cache(abspath)
        except FileNotFoundError:
            rv = path
    else:
        rv = path
    return rv


def lock_file(filepath):
    lock = filepath + ".lock"
    Path(lock).touch(exist_ok=True)


def unlock_file(filepath):
    lock = filepath + ".lock"
    if os.path.exists(lock):
        os.remove(lock)


def wait_until_unlocked_func(filepath, Nretry=10, wait=0.01):
    lock = filepath + ".lock"
    for n in range(1, Nretry+1):
        if os.path.exists(lock):
            time.sleep(wait+wait*np.random.rand())


def wait_until_unlocked(func, Nretry=10, wait=0.01):

    def inner(filepath, *args, **kwargs):
        wait_until_unlocked_func(filepath, Nretry=Nretry, wait=wait)
        return func(filepath, *args, **kwargs)

    return inner


def with_lock(func, Nretry=10, wait=0.01):

    def inner(filepath, *args, **kwargs):
        wait_until_unlocked_func(filepath, Nretry=Nretry, wait=wait)
        lock_file(filepath)
        rv = func(filepath, *args, **kwargs)
        unlock_file(filepath)
        return rv

    return inner


@wait_until_unlocked
def load_pickle(fpath, Nretry=10, wait=0.01):
    """ Load a pickle file with retries. """
    for n in range(Nretry):
        try:
            with open(fpath, "rb") as in_file:
                data = pickle.load(in_file)
            return data
        except EOFError as e:
            if n == Nretry - 1:
                raise e
        time.sleep(wait+wait*np.random.rand())


@with_lock
def save_pickle(fpath, obj):
    """ Save an object to a pickle file. """
    with open(fpath, "wb") as out_file:
        pickle.dump(obj, out_file)
