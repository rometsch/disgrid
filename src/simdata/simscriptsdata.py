# simscripts integration for data class
from .data import Data

class SData(Data):
    # data loader with simscripts support to locate remote simulations
    # and mount them via sshfs
    def __init__(self, simid, remote=True, search_args=None, **kwargs):
        self.sim = simscripts_global_lookup(simid, remote=remote, search_args=search_args)
        self.is_remote = self.sim["host"] != "localhost"
        if self.is_remote:
            self.path = self.mount()
        else:
            self.path = self.sim["path"]
        super().__init__(self.path, **kwargs)

    def mount(self):
        from simscripts.search import remote_path
        from tempfile import TemporaryDirectory
        rpath = remote_path(self.sim)
        self.tempdir = TemporaryDirectory(suffix=self.sim["uuid"], prefix="simdata-")
        lpath = self.tempdir.name
        mount_sshfs(rpath, lpath)
        return lpath

    def __del__(self):
        unmount_sshfs(self.path)
        import time
        time.sleep(0.5)

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

def mount_sshfs(remote, local):
    # mount a remote location to a local directory using sshfs
    from subprocess import run
    run(["sshfs", "-o", "ro", remote, local])

def unmount_sshfs(local):
    # unmount a sshfs mount
    from subprocess import run
    run(["fusermount", "-u", local])
