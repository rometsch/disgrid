# Integration of simdata into smurf.
# simdata is a python package to handle simulation data.
# This module provides the functionality to use sim ids
# to transparently get a data object without manually
# looking up simulation dir paths.

from smurf.search import remote_path, search
from smurf.cache import get_cache_by_id
from smurf.mount import Mount
import simdata.data


class RemoteData(simdata.data.Data):
    """ Data interface for simulations on remote hosts.
    Simdirs on remote hosts are mounted using sshfs. 
    Use user@host:path as the remote_path argument. The syntax of the remote path is equivalent to the arguments of scp."""
    def __init__(self, remote_path, cache_timeout=None, **kwargs):
        self.remote_path = remote_path
        if is_local_path(self.remote_path):
            self.path = self.remote_path
        else:
            self.path = self.mount()
        super().__init__(self.path, **kwargs)

    def mount(self, cache_timeout=None):
        self.mount_point = Mount(self.remote_path, cache_timeout=cache_timeout)
        local_path = self.mount_point.get_path()
        return local_path


class SmurfData(RemoteData):
    # data loader with smurf support to locate remote simulations
    # and mount them via sshfs
    def __init__(self, simid, search_remote=True, search_args={}, **kwargs):
        self.sim = search(simid,
                          remote=search_remote,
                          unique=True,
                          **search_args)[0]
        path = remote_path(self.sim)
        if "simdata_code" in self.sim:
            kwargs["loader"] = self.sim["simdata_code"]
        elif "simcode" in self.sim:
            kwargs["loader"] = self.sim["simcode"]
        super().__init__(path, **kwargs)
        self.simid = simid
        self.sim["simdata_code"] = self.code
        c = get_cache_by_id(self.sim["uuid"])
        c.insert(self.sim["uuid"], self.sim)


def is_local_path(path):
    """ Evaluate whether 'path' is not of the form host:path. """
    return len(path.split(":")) < 2
