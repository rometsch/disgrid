import unittest
import numpy as np
import astropy.units as u
from simdata import grid

class TestRegularGridMethods(unittest.TestCase):

    def test_passing_interfaces_only(self):
        xi1 = [1,2,3]*u.m
        xi2 = [1,2,3]*u.m
        g = grid.RegularGrid(x1_i=xi1, x2_i=xi2)
        self.assertTrue( all(g.x1_d == [1,1]*u.m))
        self.assertTrue( all(g.x1_c == [1.5, 2.5]*u.m) )
        self.assertEqual( g.dim, 2)

    def test_passing_xc_dx(self):
        xc3 = [ 1, 3.5, 7.5]*u.m
        dx3 = [ 2, 3, 5]*u.m
        g = grid.RegularGrid( x3_c = xc3, x3_d = dx3)
        self.assertTrue( all( g.x3_i == [ 0, 2, 5, 10 ]*u.m) )

    def test_no_unit(self):
        with self.assertRaises(TypeError):
            grid.RegularGrid( x1_i = np.array([1,2,3]))

    def test_fail_on_interfaces_or_centers(self):
        with self.assertRaises(ValueError):
            grid.RegularGrid( x1_d = [1,2,3]*u.m)

    def test_fail_on_only_centers(self):
        with self.assertRaises(ValueError):
            grid.RegularGrid( x1_c = [1,2,3]*u.m)


class TestCartesianGrid(unittest.TestCase):

    def test_passing_interfaces_only(self):
        xi1 = [1,2,3]*u.m
        xi2 = [1,2,3]*u.m
        g = grid.CartesianGrid(x_i=xi1, y_i=xi2)
        self.assertTrue( all(g.x_d == [1,1]*u.m))
        self.assertTrue( all(g.x_c == [1.5, 2.5]*u.m) )
        self.assertEqual( g.dim, 2)

    def test_passing_xc_dx(self):
        xc3 = [ 1, 3.5, 7.5]*u.m
        dx3 = [ 2, 3, 5]*u.m
        g = grid.CartesianGrid( z_c = xc3, z_d = dx3)
        self.assertTrue( all( g.z_i == [ 0, 2, 5, 10 ]*u.m) )


class TestCylindricalGrid(unittest.TestCase):

    def test_passing_interfaces_only(self):
        xi1 = [1,2,3]*u.m
        xi2 = [1,2,3]*u.m
        g = grid.CylindricalGrid(r_i=xi1, z_i=xi2)
        self.assertTrue( all(g.r_d == [1,1]*u.m))
        self.assertTrue( all(g.r_c == [1.5, 2.5]*u.m) )
        self.assertEqual( g.dim, 2)

    def test_passing_xc_dx(self):
        xc3 = [ 1, 3.5, 7.5]*u.m
        dx3 = [ 2, 3, 5]*u.m
        g = grid.CylindricalGrid( z_c = xc3, z_d = dx3)
        self.assertTrue( all( g.z_i == [ 0, 2, 5, 10 ]*u.m) )

class TestSphericalGrid(unittest.TestCase):

    def test_passing_interfaces_only(self):
        xi1 = [1,2,3]*u.m
        xi2 = [1,2,3]*u.m
        g = grid.SphericalGrid(r_i=xi1, theta_i=xi2)
        self.assertTrue( all(g.theta_d == [1,1]*u.m))
        self.assertTrue( all(g.theta_c == [1.5, 2.5]*u.m) )
        self.assertEqual( g.dim, 2)

    def test_passing_xc_dx(self):
        xc3 = [ 1, 3.5, 7.5]*u.m
        dx3 = [ 2, 3, 5]*u.m
        g = grid.SphericalGrid( theta_c = xc3, theta_d = dx3)
        self.assertTrue( all( g.theta_i == [ 0, 2, 5, 10 ]*u.m) )

class TestPolarGrid(unittest.TestCase):

    def test_passing_interfaces_only(self):
        xi1 = [1,2,3]*u.m
        xi2 = [1,2,3]*u.m
        g = grid.PolarGrid(r_i=xi1, phi_i=xi2)
        self.assertTrue( all(g.phi_d == [1,1]*u.m))
        self.assertTrue( all(g.phi_c == [1.5, 2.5]*u.m) )
        self.assertEqual( g.dim, 2)

    def test_passing_xc_dx(self):
        xc3 = [ 1, 3.5, 7.5]*u.m
        dx3 = [ 2, 3, 5]*u.m
        g = grid.PolarGrid( r_c = xc3, r_d = dx3)
        self.assertTrue( all( g.r_i == [ 0, 2, 5, 10 ]*u.m) )


class TestPolarGrid(unittest.TestCase):

    def test_passing_interfaces_only(self):
        xi1 = [1,2,3]*u.m
        xi2 = [1,2,3]*u.m
        g = grid.PolarGrid(r_i=xi1, phi_i=xi2)
        self.assertTrue( all(g.phi_d == [1,1]*u.m))
        self.assertTrue( all(g.phi_c == [1.5, 2.5]*u.m) )
        self.assertEqual( g.dim, 2)

    def test_passing_xc_dx(self):
        xc3 = [ 1, 3.5, 7.5]*u.m
        dx3 = [ 2, 3, 5]*u.m
        g = grid.PolarGrid( r_c = xc3, r_d = dx3)
        self.assertTrue( all( g.r_i == [ 0, 2, 5, 10 ]*u.m) )

class TestGridGetter(unittest.TestCase):

    def test_active_interfaces(self):
        xi1 = [1,2,3]*u.m
        xi2 = [1,2,3]*u.m
        g = grid.PolarGrid(r_i=xi1, phi_i=xi2, active_interfaces = ["phi"])
        self.assertTrue( all(g.get_coordinates("r") == [1.5,2.5]*u.m))
        self.assertTrue( all(g.get_coordinates("phi") == xi2) )
        self.assertEqual( g.dim, 2)

    def test_active_interfaces_wrong_definition(self):
        xi1 = [1,2,3]*u.m
        xi2 = [1,2,3]*u.m
        with self.assertRaises(ValueError):
            g = grid.PolarGrid(r_i=xi1, phi_i=xi2, active_interfaces = ["foodim"])


    def test_grid_getter(self):
        x1i = [ 1, 3.5, 7.5]*u.m
        x3c = [ 1, 3.5, 7.5]*u.m
        x3d = [ 2, 3, 5]*u.m
        g = grid.SphericalGrid( r_i = x1i, phi_c = x3c, phi_d = x3d)
        self.assertTrue( all( g.get_interfaces("r") == x1i  ))
        self.assertTrue( all( g.get_centers("phi") == x3c  ))
        self.assertTrue( all( g.get_sizes("phi") == x3d ))


if __name__ == '__main__':
    unittest.main()
