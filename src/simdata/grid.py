import numpy as np
import astropy.units as u

class Grid:
    def __init__(self, active_interfaces = [], names = []):
        self.active_interfaces = active_interfaces
        self.names = names
        for key in active_interfaces:
            if key not in names:
                raise ValueError("Can't declare active interface for dim '{}' which is not in '{}'".format(key, names))

    def get_coordinates(self, key):
        if not key in self.names:
            raise KeyError("Dimension '{}' not defined for this grid.".format(key))
        if key in self.active_interfaces:
            return self.__dict__[key+"_i"]
        else:
            return self.__dict__[key+"_c"]

    def get(self, key, coordinate_type):
        if not key in self.names:
            raise KeyError("Dimension '{}' not defined for this grid.".format(key))
        return self.__dict__[key+"_{}".format(coordinate_type)]

    def get_centers(self, key):
        return self.get(key, "c")

    def get_interfaces(self, key):
        return self.get(key, "i")

    def get_sizes(self, key):
        return self.get(key, "d")


class RegularGrid(Grid):
    def __init__(self, x1_c=None, x2_c=None, x3_c=None,
                x1_d=None, x2_d=None, x3_d=None,
                 x1_i=None, x2_i=None, x3_i=None,
                 names = ["x1", "x2", "x3"], **kwargs):
        super().__init__(**kwargs, names=names)
        # ensure units
        for n in range(1,4):
            for v in ["c", "i", "d"]:
                var = locals()["x{}_{}".format(n, v)]
                if not (isinstance(var, u.Quantity) or var is None):
                    raise TypeError("{} (type={}) doesn't have a unit".format(v+str(n), type(var)))
        count_dim = 0

        for n in range(1,4):
            xi = locals()["x{}_i".format(n)]
            xc = locals()["x{}_c".format(n)]
            dx = locals()["x{}_d".format(n)]
            if all( [xi is None, xc is None, dx is None] ):
                continue
            elif all( [xi is None, xc is None] ):
                raise ValueError("Only cell sizes given for dim nr. {}".format(n))
            elif xi is not None:
                if xc is None:
                    # calculate cell centers from cell interfaces
                    xc = 0.5*(xi[1:] + xi[:-1])
                if dx is None:
                    # calculate cell sizes from cell interfaces
                    dx = xi[1:] - xi[:-1]
            elif xc is not None:
                if xi is None and not dx is None:
                    # calculate cell interfaces from cell centers and cell sizes
                    xi = xc - 0.5*dx
                    xi = np.resize(xi, len(xi)+1).value*xi.unit
                    xi[-1] = xc[-1] + 0.5*dx[-1]
                else:
                    raise ValueError("Can't calculate cell interfaces for dim nr. {}. Cell sizes dx not given".format(n))


            count_dim += 1
            self.__dict__["{}_c".format(names[n-1])] = xc
            self.__dict__["{}_i".format(names[n-1])] = xi
            self.__dict__["{}_d".format(names[n-1])] = dx

        self.dim = count_dim

def coordinate_args(loc, names):
    types = ["c", "i", "d"]
    old = ["{}_{}".format(a,b) for a in names for b in types]
    new = ["{}_{}".format(a,b) for a in ["x1", "x2", "x3"] for b in types]
    rv = {}
    for o, n in zip(old, new):
        rv[n] = loc[o]
    return rv


class CartesianGrid(RegularGrid):
    def __init__(self, x_c=None, y_c=None, z_c=None,
                 x_i=None, y_i=None, z_i=None,
                 x_d = None, y_d=None, z_d=None, **kwargs):
        names = ["x", "y", "z"]
        super().__init__(**coordinate_args(locals(), names), **kwargs, names = names)

class CylindricalGrid(RegularGrid):
    def __init__(self, r_c=None, phi_c=None, z_c=None,
                 r_i=None, phi_i=None, z_i=None,
                 r_d = None, phi_d=None, z_d=None, **kwargs):
        names = ["r", "phi", "z"]
        super().__init__(**coordinate_args(locals(), names), **kwargs, names = names)

class SphericalGrid(RegularGrid):
    def __init__(self, r_c=None, phi_c=None, theta_c=None,
                 r_i=None, phi_i=None, theta_i=None,
                 r_d = None, phi_d=None, theta_d=None, **kwargs):
        names = ["r", "phi", "theta"]
        super().__init__(**coordinate_args(locals(), names), **kwargs, names = names)

class PolarGrid(RegularGrid):
    def __init__(self, r_c=None, phi_c=None,
                 r_i=None, phi_i=None,
                 r_d = None, phi_d=None, **kwargs):
        names = ["r", "phi"]
        super().__init__(**coordinate_args(locals(), names), **kwargs, names = names)

