import astropy.units as u
from .grid import Grid

class Field:
    def __init__(self, grid, data, time, name):
        for x, n in zip((data, time), ["data", "time"]):
            ensure_is_astropy_quantity(x, name=n)

        if not isinstance(grid, Grid):
            raise TypeError("grid is not a valid Grid class (type={})".format(type(grid)))

        self.grid = grid
        self.data = data
        self.time = time
        self.name = name

    def get_grid(self):
        return self.grid

    def get_data(self):
        return self.data

    def get_name(self):
        return self.name

    def get_time(self):
        return self.time

    def get(self, key):
        return self.__dict__[key]

def ensure_is_astropy_quantity(x, name="variable"):
    # ensure data and time have astropy unit=
    if not (isinstance(x, u.Quantity) or isinstance(x, u.Unit)):
        raise TypeError("{} (type={}) doesn't have a unit".format(n, type(x)))
