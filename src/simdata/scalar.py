# scalar quantities and time series
# data is stored as Ntimes x Ndata
from .time_access import get_indices_time_interval, get_index_closest_time, ensure_unit


class Scalar:
    def __init__(self, time, data, name=""):
        self.time = time
        self.data = data
        self.name = name
        for v in [time, data]:
            ensure_unit(v)

    def __getitem__(self, key):
        return self.get(key)

    def get(
        self,
        key=None,
        t=None,
        tmin=None,
        tmax=None,
        return_time=False,
    ):
        """ Extract a range of data.

        How data is extracted depends on which keyword is set.
        Using the key parameter, the data and time arrays are simply indexed.
        When using the time parameters, the closest matching point in time is
        caculated and data is returned for this closest matching time.
        When using tmin and tmax parameters, all datapoints belonging to points
        in time (t) with tmin <= t <= tmax are extracted.

        Parameters
        ----------
        key : integer or list like or range or slice
            Key which which an array can be accessed.
        t : :obj:`astropy.units.quantity.Quantity`
            Time-like quantity for which the closest matching value is selected.
        tmin : :obj:`astropy.units.Quantity`
            Lower bound on time.
        tmax : :obj:`astropy.units.Quantity`
            Upper bound on time.
        return_time : bool
            Whether or not to return a time array alongside the data.

        Returns
        -------
        :obj:`astropy.units.Quantity` or tuple of such
            Return the data or if return_time is true a tuple with (time, data).
        """
        rv = self.data
        time = self.time
        if key is not None:
            rv = rv[key]
            time = time[key]
        if tmin is not None and tmax is not None:
            inds = get_indices_time_interval(tmin, tmax, time)
            rv = rv[inds]
            time = time[inds]
        if t is not None:
            ind = get_index_closest_time(t, time)
            rv = rv[ind]
            time = time[ind]
        if return_time:
            rv = (time, rv)
        return rv

    def get_time_interval(self, tmin, tmax, **kwargs):
        inds = get_indices_time_interval(tmin, tmax, self.time)
        return self.get(inds, **kwargs)

    def get_closest_to_time(self, t, **kwargs):
        """ Get the data that's closest to the given time.

        Parameters
        ----------
        t: :obj:`astropy.units.quantity.Quantity`
            Time for which to get the value.

        Returns:
        :obj:`astropy.units.Quantity`
            Datapoint for which time is closes to t.
        """
        ind = get_index_closest_time(t, self.time)
        return self.get(key=ind, **kwargs)
