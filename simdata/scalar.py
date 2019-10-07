# scalar quantities and time series
# data is stored as Ntimes x Ndata
import numpy as np
import astropy.units as u

class Scalar:

    def __init__(self, time, data, axes=[], time_dim=0, axes_dim=None, name=""):
        self.time = time
        self.data = data
        self.name = name
        self.axes = axes
        self.axes_dim = axes_dim
        self.time_dim = time_dim
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
        direct_key = False
        if key is None:
            if low is None and up is None:
                inds = slice(0, self.data.shape[self.time_dim])
            elif low is not None and up is not None:
                inds = slice(low,up)
            elif low is None:
                inds = slice(0, up)
            elif up is None:
                inds = slice(low, self.data.shape[self.time_dim])
        elif type(key) in [range, slice]:
            inds = key
        elif type(key) in [int, np.int64, np.int32, np.int16]:
            inds = int(key)
        else:
            try:
                len(key)
                inds = key
                direct_key = True
            except TypeError:
                try:
                    inds = int(key)
                    direct_key = True
                except TypeError:
                    raise TypeError("Key '{}' (type='{}') is neither a range nor a list of indices".format(key, type(key)))
        rv = None
        if axis is not None:
            # return a specific axis
            if axis in self.axes:
                n_axis = (n for n,s in enumerate(self.axes) if s == axis).__next__()
            else:
                raise KeyError("Axis '{}' not in available axes '{}'".format(axis, self.axes))
        else:
            n_axis = slice(0, len(self.axes))
        # select the correct values
        if direct_key:
            rv = self.data[key]
        elif self.time_dim == 0:
            if self.axes_dim is None:
                rv = self.data[inds]
            else:
                rv = self.data[inds, n_axis]
        elif self.time_dim == 1:
            if self.axes_dim is None:
                rv = self.data[:,inds]
            else:
                rv = self.data[:, inds, n_axis]
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

    def get_time_closest(self, t, **kwargs):
        if not isinstance(t, u.Quantity):
            raise TypeError("'{}' does not have a unit!".format(t))
        ind = np.argmin( np.abs(self.time-t))
        return self.get(key=ind, **kwargs)
