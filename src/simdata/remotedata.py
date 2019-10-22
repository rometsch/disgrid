# sshfs integration for data class
from .data import Data
from .mount import Mount

class RemoteData(Data):
    """ Data interface for simulations on remote hosts.
    Simdirs on remote hosts are mounted using sshfs. 
    Use user@host:path as the remote_path argument. The syntax of the remote path is equivalent to the arguments of scp."""
    def __init__(self, remote_path, **kwargs):
        self.remote_path = remote_path
        if is_local_path(self.remote_path):
            self.path = self.remote_path
        else:
            self.path = self.mount()
        super().__init__(self.path, **kwargs)

    def mount(self):
        self.mount_point = Mount(self.remote_path)
        local_path = self.mount_point.get_path()
        return local_path

def is_local_path(path):
    """ Evaluate whether 'path' is not of the form host:path. """
    return len(path.split(":")) < 2
