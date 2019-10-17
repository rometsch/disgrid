# simscripts integration for data class
from .data import Data
from simscripts.search import remote_path
from .mount import Mount

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
        self.mount = Mount(remote_path(self.sim), self.sim["uuid"])
        local_path = self.mount.get_path()
        return local_path

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

