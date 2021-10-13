# Integration of the simdata-net package.
# This plugin adds the NData class which uses the server-client model to get data.
# Online the get method is supported.

import os
import urllib
from simdata_net.client import simdata_request
import diskcache
from simdata.config import Config

class NData:
    
    def __init__(self, simid, caching=True):
        self.simid = simid
        self.caching = caching
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
            size_limit = 10*8**12 # 10 GiB
            self.config["cachesize"] = size_limit
            
        self.cache = diskcache.Cache(directory=cachedir, 
                                     eviction_policy="least-recently-used",
                                     size_limit=size_limit)
    
    def get(self, var=None, dim=None, N=None, planet=None, t=None, fluid="gas", **kwargs):
        """ General getter function. """
        url = urllib.parse.urlencode(dict(
            simid=self.simid,
            var=var, dim=dim, N=N, planet=planet, t=t, fluid=fluid
        ))
        try:
            rv = self.cache[url]
        except KeyError:
            rv = simdata_request(url)
            if self.caching:
                self.cache.add(url, rv)
        return rv