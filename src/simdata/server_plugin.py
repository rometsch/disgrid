# Integration of the simdata-net package.
# This plugin adds the NData class which uses the server-client model to get data.
# Online the get method is supported.

from .config import Config

from simdata_net.client import simdata_request


class NData:

    def __init__(self, simid):
        self.simid = simid
        self.config = Config()
        
        try:
            self.relay = self.config.data["relay-server"]
        except KeyError:
            self.relay = None

    def get(self, var=None, dim=None, N=None, planet=None, t=None, fluid="gas", **kwargs):
        """ General getter function. """
        query = dict(var=var, dim=dim, N=N, planet=planet,
                     t=t, fluid=fluid, **kwargs)

        rv = simdata_request(self.simid, query, hostname=self.relay)
        return rv
