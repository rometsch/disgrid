# Integration of simdata into smurf.
# simdata is a python package to handle simulation data.
# This module provides the functionality to use sim ids
# to transparently get a data object without manually
# looking up simulation dir paths.

import os
import pickle

from . import data
from .config import Config
from smurf.cache import LocalSimCache, get_cache_by_id, CacheMiss
from smurf.info import Info
from smurf.mount import Mount
from smurf.search import remote_path, search


class RemoteData(data.Data):
    """ Data interface for simulations on remote hosts.
    Simdirs on remote hosts are mounted using sshfs. 
    Use user@host:path as the remote_path argument. The syntax of the remote path is equivalent to the arguments of scp."""

    def __init__(self, remote_path, cache_timeout=None, **kwargs):
        self.remote_path = remote_path
        if is_local_path(self.remote_path):
            self.path = self.remote_path
            path = self.path
        # else:
        #     self.path = self.mount()
        else:
            path = None
        super().__init__(path, **kwargs)

    def init(self):
        if self.path is None:
            self.mount()
        super().init()

    def mount(self, cache_timeout=None):
        remote_path_root = self.remote_path.split(":")[0] + ":/"
        remote_path_child = self.remote_path.split(":")[1].lstrip("/")
        self.mount_point = Mount(remote_path_root, cache_timeout=cache_timeout)
        local_path = os.path.join(
            self.mount_point.get_path(), remote_path_child)
        self.path = local_path


class SmurfData(RemoteData):
    # data loader with smurf support to locate remote simulations
    # and mount them via sshfs
    def __init__(self, simid, search_remote=True, search_args={}, **kwargs):
        simid = insert_local_sim_to_cache(simid)
        self.sim = search(simid,
                          remote=search_remote,
                          unique=True,
                          **search_args)[0]
        path = remote_path(self.sim)
        if "simdata_code" in self.sim:
            kwargs["loader"] = self.sim["simdata_code"]
        elif "simcode" in self.sim:
            kwargs["loader"] = self.sim["simcode"]
        self.simid = simid
        self.config = Config()
        spec = self.search_spec()
        if spec is not None:
            kwargs["spec"] = spec
        super().__init__(path, init_hooks=[self.register_code, self.save_spec], **kwargs)

    def register_code(self):
        self.sim["simdata_code"] = self.code
        c = get_cache_by_id(self.sim["uuid"])
        c.insert(self.sim["uuid"], self.sim)

    def search_spec(self):
        spec = None
        simid = self.sim["uuid"]
        if "specdir" in self.config.data:
            specdir = self.config.data["specdir"]
        elif "cachedir" in self.config.data["cachedir"]:
            specdir = self.config.data["cachedir"]
        else:
            spec = None
        specfile = os.path.join(specdir, f"{simid}.spec.pickle")
        if os.path.exists(specfile):
            try:
                with open(specfile, "rb") as in_file:
                    spec = pickle.load(in_file)
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
        print(self.config.data)
        if "specdir" in self.config.data:
            specdir = self.config.data["specdir"]
        elif "cachedir" in self.config.data["cachedir"]:
            specdir = self.config.data["cachedir"]
        else:
            return
        if hasattr(self.loader, "spec"):
            specfile = os.path.join(specdir, f"{simid}.spec.pickle")
            with open(specfile, "wb") as out_file:
                pickle.dump(self.loader.spec, out_file)


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


def is_local_path(path):
    """ Evaluate whether 'path' is not of the form host:path. """
    return len(path.split(":")) < 2
