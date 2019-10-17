from subprocess import run, PIPE
import sys
import time
import os
import atexit
from tempfile import mkdtemp
from pathlib import Path
from .doublefork import detachify
from uuid import uuid4

class Mount:
    def __init__(self, remote, uuid):
        self.remote = remote
        self.sim_uuid = uuid
        self.mount_uuid = str(uuid4())
        self.mount()

    def mount(self):
        existing = find_existing_mount(self.remote)
        if existing:
            self.reuse_mount(existing.parent)
        else:
            self.create_mount()
            self.create_unmounter()
        self.flag_active()
        self.flag_finished_on_exit()

    def reuse_mount(self, path):
        self.tempdir = path
        self.flag_active()

    def create_mount(self):
        self.tempdir = Path( mkdtemp(suffix="-"+self.sim_uuid, prefix="simdata-") )
        os.mkdir( self.tempdir / "active" )
        os.mkdir( self.get_path() )
        mount_sshfs(self.remote, self.get_path(), remove=True)

    def flag_active(self):
        """ Flag the mount point to be in use. """
        self.get_active_flag_file().touch()

    def get_active_flag_file(self):
        return self.tempdir / "active" / self.mount_uuid

    def flag_finished(self):
        """ Flag the mount point to be in not in use anymore. """
        self.get_finished_flag_file().touch()
        os.remove( self.get_active_flag_file() )

    def time_since_finished(self):
        """ Return the time since the finished flag file was last touched """
        if self.get_finished_flag_file().exists():
            return time.time() - self.get_finished_flag_file().stat().st_ctime
        else:
            return -1

    def in_use(self):
        """ Check whether there are still users of the mount point """
        return len( os.listdir( self.tempdir / "active" ) ) > 0

    def get_finished_flag_file(self):
        return Path( self.tempdir ) / "finished"

    def flag_finished_on_exit(self):
        atexit.register( lambda : self.flag_finished() )

    def create_unmounter(self):
        """ Spawn a detached process to unmount the sshfs after a certain period. """
        atexit.register( lambda : unmount_delayed(self, remove=True) )

    def get_path(self):
        """ Return the path of the mount point """
        return os.path.join(self.tempdir, "mnt")

    def unmount(self):
        """ Unmount and remove temporary directories """
        unmount_sshfs( self.get_path() )
        os.rmdir( self.get_path() )
        os.rmdir( self.tempdir / "active" )
        os.remove( self.get_finished_flag_file() )
        os.rmdir( self.tempdir )

def mount_sshfs(remote, local, remove=True):
    # mount a remote location to a local directory using sshfs
    run(["sshfs", "-o", "ro", remote, local])

def find_existing_mount(path):
    res = run(["mount"], stdout=PIPE, encoding="utf-8").stdout.splitlines()
    for line in res:
        if path in line:
            return Path(line.split()[2])

def unmount_sshfs(local):
    # unmount a sshfs mount
    if sys.platform == "darwin":
        run(["umount", "-f", local])
    else:
        run(["fusermount", "-u", local])

@detachify
def unmount_delayed(mount, remove=False):
    """ unmount after some time if not still in use """
    time.sleep(1)
    waiting_time = 120
    # check whether path is still mounted
    for n in range(1000):
        if not mount.in_use() and mount.time_since_finished() > waiting_time:
            mount.unmount()
            return
        else:
            time.sleep(waiting_time)
