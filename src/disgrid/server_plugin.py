# Integration of the disgrid-net package.
# This plugin adds the NData class which uses the server-client model to get data.
# Online the get method is supported.

import urllib
from smurfnet.client import make_request
from disgrid.config import Config


class NData:

    def __init__(self, simid, update=False):
        self.simid = simid
        self.config = Config()
        self.update = update


    def get(self, update=False, var=None, dim=None, N=None, planet=None, t=None, fluid=None, **kwargs):
        """ General getter function. """
        query = dict(
            simid=self.simid,
            var=var, dim=dim, N=N, planet=planet, t=t, fluid=fluid
        )
        query = {k: v for k, v in query.items() if v is not None}
        uri = urllib.parse.urlencode(query)
        url = f"disgrid://localhost/get?{uri}"
        if update or self.update:
            url += "#update"
        rv = make_request(url)
        return rv

    def avail(self, update=False):
        """ Get the available fields. """
        query = dict(simid=self.simid)
        uri = urllib.parse.urlencode(query)
        url = f"disgrid://localhost/avail?{uri}"
        if update or self.update:
            url += "#update"
        rv = make_request(url)
        return rv
