


class Grid:
    def __init__(self, dim):
        self.dim = dim


class CartesianGrid(Grid):
    def __init__(self, dim, x=None, y=None, z=None,
                dx=None, dy=None, dz=None, 
                xi=None, yi=None, zi=None):
        if dim != sum([val is not None for val in [x, y, z]]):
            raise ValueError("Dimension {} does not fit passed arguments {} {} {}".format(dim, x, y, z))
        self.dim = dim



class CartesianGrid(Grid):
    def __init__(self, dim)
