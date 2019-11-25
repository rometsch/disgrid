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
    ensure_in_time_interval(t, time)
    ind = np.argmin(np.abs(time - t))
    return ind


def ensure_unit(x):
    """ Ensure that x is an astropy quantity. Raise an error otherwise. """
    if not isinstance(x, u.Quantity):
        raise TypeError("'{}' does not have a unit!".format(x))


def ensure_in_time_interval(t, time, delta=0.001):
    """ Ensure that t is in the time interval of time.

    Raise an error if time is not within the interval defined by
    time +- some tolerance defined by delta.

    Parameters
    ----------
    t : :obj:`astropy.units.quantity.Quantity`
        Time to be checked.
    time : :obj:`astropy.units.quantity.Quantity`
        Time array to define the allowed interval.
    delta : float
        Percentage of the time interval to add as tolerance.

    Raises
    ------
    ValueError
        If t is not inside time interval + tolerance.
    """
    tmin = time[0] - delta * (time[-1] - time[0])
    tmax = time[-1] + delta * (time[-1] - time[0])
    if t > tmax or t < tmin:
        raise ValueError(
            "t={} is not inside time interval [{}, {}]".format(t, tmin, tmax),
            ", delta = {}".format(delta))
