# vector quantities and time series
# data is stored as Ntimes x Ndata
import numpy as np
import astropy.units as u

class Vector:

    def __init__(self, time, data, axes=[], name=""):
        self.time = time
        self.data = data
        self.name = name
        self.axes = axes
        for v in [time, data]:
            if not isinstance(v, u.Quantity):
                raise TypeError("'{}' is a physical value but not an astropy quantity!".format(v))

    def __getitem__(self, key):
        return self.get(key)

    def get(self, key=None, low=None, up=None, axis=None, return_time=False):
        # TODO: handle keys
        # handle
        # 1) a single number
        # 2) array of indices [1, 5, 7 ,8]
        # 3) range operator range(2,20)
        if key is None:
            if low is None and up is None:
                raise ValueError("No key and neither upper nor lower bound given!")
            if low is not None and up is not None:
                inds = range(low,up)
            elif low is None:
                inds = range(0, up)
            elif up is None:
                inds = range(low, len(self.data))
        elif type(key) in [int, range, slice]:
            inds = key
        else:
            try:
                len(key)
                inds = key
            except TypeError:
                raise TypeError("Key is neither a range nor a list of indices")
        rv = None
        if axis is not None:
            if axis in self.axes:
                n_axis = (n for n,s in enumerate(self.axes) if s == axis).__next__()
                rv = self.data[inds, n_axis]
        else:
            rv = self.data[inds]
        if return_time:
            rv = (self.time[inds], rv)
        return rv

    def get_time_interval(self, tmin, tmax, **kwargs):
        for t in [tmin, tmax]:
            if not isinstance(t, u.Quantity):
                raise TypeError("'{}' does not have a unit!".format(t))
        inds = np.logical_and( tmin <= self.time, tmax >= self.time)
        key = np.arange(len(inds))[inds]
        return self.get(key, **kwargs)
