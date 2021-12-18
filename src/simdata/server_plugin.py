# Integration of the simdata-net package.
# This plugin adds the NData class which uses the server-client model to get data.
# Online the get method is supported.

import os
import urllib
from smurfnet.client import make_request
import diskcache
from simdata.config import Config


class NData:

    def __init__(self, simid, caching=True, update=False):
        self.simid = simid
        self.caching = caching
        self.config = Config()
        self.update = update

        try:
            self.relay = self.config.data["relay-server"]
        except KeyError:
            self.relay = "localhost"

        self._init_cache_()

    def _init_cache_(self):
        if not self.caching:
            self.cache = dict({})
            return

        self.config = Config()
        try:
            cachedir = self.config["cachedir"]
        except KeyError:
            cachedir = os.path.expanduser("~/.simdata/cache")
            self.config["cachedir"] = cachedir

        try:
            size_limit = self.config["cachesize"]
        except KeyError:
            size_limit = 10*8**12  # 10 GiB
            self.config["cachesize"] = size_limit

        self.cache = diskcache.Cache(directory=cachedir,
                                     eviction_policy="least-recently-used",
                                     size_limit=size_limit)

    def get(self, update=False, var=None, dim=None, N=None, planet=None, t=None, fluid=None, **kwargs):
        """ General getter function. """
        query = dict(
            simid=self.simid,
            var=var, dim=dim, N=N, planet=planet, t=t, fluid=fluid
        )
        query = {k: v for k, v in query.items() if v is not None}
        uri = urllib.parse.urlencode(query)
        try:
            if update or self.update:
                raise KeyError
            rv = self.cache[f"get?{uri}"]
        except KeyError:
            url = f"simdata://{self.relay}/get?{uri}"
            if update or self.update:
                url += "#update"
            rv = make_request(url)
            if self.caching:
                self.cache.add(f"get?{uri}", rv)
        return rv

    def avail(self, update=False):
        """ Get the available fields. """
        query = dict(simid=self.simid)
        uri = urllib.parse.urlencode(query)
        try:
            if update or self.update:
                raise KeyError
            rv = self.cache[f"avail?{uri}"]
        except KeyError:
            url = f"simdata://{self.relay}/avail?{uri}"
            if update or self.update:
                url += "#update"
            rv = make_request(url)
            if self.caching:
                self.cache.add(uri, rv)
        return rv
