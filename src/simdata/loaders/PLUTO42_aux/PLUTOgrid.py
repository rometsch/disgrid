import numpy as np


class Grid:
    """grid structure
    features:
            DIR:    direction, can be (0, 1, 2) (r, phi, blank for 2D polar; r, theta, phi for 3D spherical, etc)
            Ncells: number of cells in given direction
            x:      cell-centered values
            xi:     cell interfaces
            dx:     cell size
            xl:     left interfaces
            xr:     right interfaces

    Input
        DIR: the direction in which the grid is constructed (e.g. 0 for 'r' in a polar setup, 1 for 'y' in a cartesian setup)
        origin: the directory where 'grid.out' is searched for
    Output
        None
    
    Example
        gr = Grid(0, origin='/some/data/folder/containing/gridout/')
        r = gr.x

        g = [Grid(i, origin='/some/data/folder/containing/gridout/') for i in [0,1,2]]
        dphi = g[1].dx

    credit goes to sdoetsch for inspiration on improving the file-reading part: https://gitlab.mpcdf.mpg.de/sdoetsch/plutoplot.git
    """
    def __init__(self, DIR=0, origin='.'):
        filename = "%s/grid.out" % origin
        try:
            file = open(filename)
            file.close()
        except:
            raise FileNotFoundError(
                'PLUTO dependencies: No such file: \'%s\'' % filename)

        dim_counter = -1
        with open(filename, 'r') as f:
            while True:
                line = f.readline()
                if line[0] == '#': continue  #ignore comments
                splitted = line.split()
                if len(splitted) == 1:
                    dim_counter += 1
                    NXi = int(splitted[0])
                    data = np.fromfile(f, sep=' ',
                                       count=NXi * 3).reshape(NXi, 3)
                    if dim_counter != DIR: continue
                    self.xl = data[:, 1]
                    self.xr = data[:, 2]
                    self.x = 0.5 * (self.xl + self.xr)
                    self.xi = np.append(self.xl, self.xr[-1])
                    self.dx = self.xr - self.xl
                    self.Ncells = NXi
                    self.DIR = DIR
                    return


def resolve_geometry(origin='.'):
    """looks inside a directory containing PLUTO data (specifically a grid.out file)
    and extracts information about the geometry of the problem
        Input
            origin: the folder to look inside
        Output
            dimensions: number of dimensions of the problem (1, 2, 3).
            geometry:   can be "CARTESIAN", "POLAR", "CYLINDRICAL", "SPHERICAL".
            coords:     PLUTO is a bit ambiguous about its coord system at first glance. Detailed info on the coord system is included here.
    """
    dimensions, geometry = -1, "NONE"  #defaults
    with open(origin + '/grid.out', 'r') as f:
        while True:
            l = f.readline()
            if l.startswith('# DIMENSIONS'):
                dimensions = int(l.strip().split()[-1])
            if l.startswith('# GEOMETRY'):
                geometry = l.strip().split()[-1]
                break
    if dimensions == -1 or geometry == "NONE":
        raise DataNotFoundError("Could not resolve geometry from grid.out.")

    if dimensions == 1:
        if geometry == "CARTESIAN": coords = "x"
        elif geometry in ["POLAR", "CYLINDRICAL", "SPHERICAL"]: coords = "rad"
        else:
            raise NotImplementedError(
                "I don't understand this geometry: DIMENSIONS = %d, GEOMETRY = %s."
                % (dimensions, geometry))
    if dimensions == 2:
        if geometry == "CARTESIAN": coords = "x y"
        elif geometry == "POLAR": coords = "rad azimuth"
        elif geometry == "CYLINDRICAL": coords = "rad z"
        elif geometry == "SPHERICAL": coords = "rad polar"
        else:
            raise NotImplementedError(
                "I don't understand this geometry: DIMENSIONS = %d, GEOMETRY = %s."
                % (dimensions, geometry))
    elif dimensions == 3:
        if geometry == "CARTESIAN": coords = "x y z"
        elif geometry == "POLAR": coords = "rad azimuth z"
        elif geometry == "SPHERICAL": coords = "rad polar azimuth"
        else:
            raise NotImplementedError(
                "I don't understand this geometry: DIMENSIONS = %d, GEOMETRY = %s."
                % (dimensions, geometry))

    return [dimensions, geometry, coords.split()]
