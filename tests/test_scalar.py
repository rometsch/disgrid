import unittest
import astropy.units as u
import numpy as np
from simdata import scalar


class TestDataMethods(unittest.TestCase):
    def setUp(self):
        self.t = 0.5 * np.arange(10) * u.s
        np.random.seed(0)
        self.d = np.random.rand(10) * u.cm
        name = 'distance'
        self.v = scalar.Scalar(self.t, self.d, name=name)

    def test_no_unit(self):
        with self.assertRaises(TypeError):
            scalar.Scalar([0, 1], [1, 2] * u.s)
        with self.assertRaises(TypeError):
            scalar.Scalar([0, 1] * u.s, [1, 2])

    def test_get_ind_list(self):
        self.assertTrue(all(self.v.get([2, 4, 7]) == self.d[[2, 4, 7]]))

    def test_get_range(self):
        self.assertTrue(all(self.v.get(range(2, 8)) == self.d[2:8]))

    def test_get_slice(self):
        self.assertTrue(all(self.v.get(slice(2, 8, 2)) == self.d[2:8:2]))
        self.assertTrue(all(self.v[2:8:2] == self.d[2:8:2]))

    def test_get_time(self):
        t, v = self.v.get(range(2, 8), return_time=True)
        self.assertTrue(all(t == self.t[2:8]))

    def test_get_time_interval(self):
        t, v = self.v.get_time_interval(0.75 * u.s,
                                        3.9 * u.s,
                                        return_time=True)
        self.assertTrue(all(t == self.t[2:8]))
        with self.assertRaises(TypeError):
            self.v.get_time_interval(0.75, 3.9 * u.s, return_time=True)

    def test_get_time_interval_via_get(self):
        t, v = self.v.get(tmin=0.75 * u.s, tmax=3.9 * u.s, return_time=True)
        self.assertTrue(all(t == self.t[2:8]))
        with self.assertRaises(TypeError):
            self.v.get(tmin=0.75, tmax=3.9 * u.s, return_time=True)

    def test_get_time_closest_via_get(self):
        t, v = self.v.get(t=0.75 * u.s, return_time=True)
        self.assertTrue(t == self.t[1])
        with self.assertRaises(TypeError):
            self.v.get(t=0.75, return_time=True)

    def test_get_by_time_out_of_interval(self):
        t, v = self.v.get(t=4.5001 * u.s, return_time=True)
        self.assertTrue(t == self.t[-1])
        with self.assertRaises(ValueError):
            self.v.get(t=5*u.s, return_time=True)

if __name__ == '__main__':
    unittest.main()
