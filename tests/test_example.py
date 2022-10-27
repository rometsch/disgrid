import os
import unittest
from disgrid import data
from disgrid.loaders import example_loader
import astropy.units as u
import numpy as np

from run import get_repo_abspath
code_sample_path = os.path.join(get_repo_abspath(), "samples/example")


class TestExampleLoader(unittest.TestCase):
    def setUp(self):
        self.d = data.Data(code_sample_path)
        self.d.init()

    def test_identify_code_directly(self):
        self.assertTrue(example_loader.identify(code_sample_path))

    def test_identify_code_via_class(self):
        loader = example_loader.Loader(code_sample_path)
        self.assertEqual(loader.code_info, ("example", "0.0", "description"))

    def test_identify_code_via_data(self):
        self.assertEqual(self.d.loader.code_info, ("example", "0.0", "description"))

    def test_params(self):
        self.assertEqual(self.d.parameters["tfinal"], 100.0)

    def test_time(self):
        self.assertEqual(self.d.loader.get_time(3), 1500*u.s)

    def test_domain_size(self):
        self.assertEqual(self.d.loader.get_time(3), 1500*u.s)

    def test_planet(self):
        x = self.d.get("x", planet=0)
        self.assertEqual(x[3], -0.5480077554195742 * 2 * u.cm)

    def test_fluids(self):
        f = self.d.fluids["He"]
        self.assertEqual(f.name, "He")

    def test_scalar(self):
        x = self.d.get(fluid="He", var="momentum", dim="scalar")
        self.assertEqual(x[2], 0.25*(2*u.cm)/(20*u.s))

    def test_2d_shape(self):
        vx = self.d.get(fluid="He", var="vx", dim="2d", N=2)
        self.assertEqual(vx.data.shape, (11,10))

    def test_2d_grid(self):
        vx = self.d.get(fluid="He", var="vx", dim="2d", N=2)
        self.assertEqual(vx.grid.get_coordinates("x")[1], -0.8*2*u.cm)
        
    def test_2d_time(self):
        vx = self.d.get(fluid="He", var="vx", dim="2d", N=2)
        self.assertEqual(vx.time, 1000*u.s)
        


if __name__ == '__main__':
    unittest.main()
