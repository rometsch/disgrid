# simscripts integration for data class
from subprocess import run
from .data import Data
import sys
import time

class SData(Data):
    # data loader with simscripts support to locate remote simulations
    # and mount them via sshfs
    def __init__(self, simid, remote=True, search_args=None, **kwargs):
        self.sim = simscripts_global_lookup(simid, remote=remote, search_args=search_args)
        self.search_remote = remote
        self.is_remote = self.sim["host"] != "localhost"
        if self.is_remote:
            self.path = self.mount()
        else:
            self.path = self.sim["path"]
        super().__init__(self.path, **kwargs)

    def mount(self):
        from simscripts.search import remote_path
        from tempfile import mkdtemp
        rpath = remote_path(self.sim)
        self.tempdir = mkdtemp(suffix=self.sim["uuid"], prefix="simdata-")
        lpath = self.tempdir
        mount_sshfs(rpath, lpath, remove=True)
        return lpath

def simscripts_global_lookup(pattern, remote=True, search_args=None):
    try:
        from simscripts.search import search
        try:
            if search_args is None:
                search_args = {}
            rv = search(pattern, remote=remote, **search_args)
            # use the first match
            if len(rv) == 0:
                raise KeyError("Could not locate simulation for pattern : '{}'".format(pattern))
            return rv[0]

        except Exception as e:
            print("Simscripts lookup failed with exception: {}".format(e))
            raise
    except ImportError:
        print("Simscripts is not installed on this maschine!")
        return {}

def mount_sshfs(remote, local, remove=True):
    # mount a remote location to a local directory using sshfs
    run(["sshfs", "-o", "ro", remote, local])
    import atexit
    atexit.register( lambda : unmount_sshfs(local, remove=remove) )

def unmount_sshfs(local, remove=False):
    # unmount a sshfs mount
    if sys.platform == "darwin":
        run(["umount", "-f", local])
    else:
        run(["fusermount", "-u", local])
    import os
    os.rmdir(local)
