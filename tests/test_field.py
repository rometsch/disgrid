import unittest
import importlib
import astropy.units as u
from disgrid import field
from disgrid import grid


class TestFluidMethods(unittest.TestCase):
    def setUp(self):
        pass

    def test_no_unit(self):
        with self.assertRaises(TypeError):
            grid.Grid.__init__ = lambda x: None
            field.Field(grid.Grid(), 2, 3, "test")
        importlib.reload(grid)

    def test_has_unit(self):
        grid.Grid.__init__ = lambda x: None
        field.Field(grid.Grid(), 2 * u.m, 3 * u.m, "test")
        importlib.reload(grid)

    def test_not_grid_type(self):
        with self.assertRaises(TypeError):
            field.Field(1, 2 * u.m, 3 * u.m, "test")


if __name__ == '__main__':
    unittest.main()
