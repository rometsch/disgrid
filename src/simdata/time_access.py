# Helper functions to handle accessing array by time.
import astropy.units as u
import numpy as np


def get_indices_time_interval(tmin, tmax, time):
    """ Find the indeces for which points in time fullfil tmin <= t <= tmax.
    
    Parameters
    ----------
    tmin : :obj:`astropy.units.Quantity`
            Lower bound on time.
        tmax : :obj:`astropy.units.Quantity`
            Upper bound on time.
    time : :obj:`astropy.units.quantity.Quantity`
        Time array to find the index for.

    Returns
    -------
    array of integers
        Index of the closest point in time.
    """
    for t in [tmin, tmax, time]:
        ensure_unit(t)
    mask = np.logical_and(tmin <= time, tmax >= time)
    inds = np.arange(len(mask))[mask]
    return inds


def get_index_closest_time(t, time):
    """ Find the index for the point in time thats closest to t.
    
    Parameters
    ----------
    t : :obj:`astropy.units.quantity.Quantity`
        Time for which to find the index
    time : :obj:`astropy.units.quantity.Quantity`
        Time array to find the index for.

    Returns
    -------
    integer
        Index of the closest point in time.
    """
    ensure_unit(time)
    ensure_unit(t)
    ind = np.argmin(np.abs(time - t))
    return ind


def ensure_unit(x):
    """ Ensure that x is an astropy quantity. Raise an error otherwise. """
    if not isinstance(x, u.Quantity):
        raise TypeError("'{}' does not have a unit!".format(x))
