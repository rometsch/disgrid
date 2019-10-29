import unittest
from simdata import data
from simdata.loaders import PLUTO42
import astropy.units as u
import astropy.constants as const
import numpy as np
import os

from run import get_repo_abspath
code_sample_path = os.path.join(get_repo_abspath(), "samples/PLUTO42/out")


class TestPLUTO42Loader(unittest.TestCase):
    def setUp(self):
        self.d = data.Data(code_sample_path)

    def test_identify_code_via_data(self):
        self.assertEqual(self.d.loader.code_info, ("PLUTO", "4.2", "vanilla"))

    def test_identify_code_directly(self):
        self.assertTrue(PLUTO42.identify(code_sample_path))

    def test_units(self):
        self.assertTrue(
            close(1. * self.d.loader.units['length'],
                  5.2 * 1.49597871e13 * u.cm))
        self.assertTrue(close(1. * self.d.loader.units['mass'], 2e33 * u.g))
        self.assertTrue(
            close(
                1. * self.d.loader.units['time'],
                2 * np.pi * np.sqrt((5.2 * 1.49597871e13)**3 * u.cm**3 /
                                    (const.G.cgs * 2e33 * u.g)).to(u.s)))

    def test_fine_output_times(self):
        #initial time is 0 -> removed first element
        time = (self.d.fluids["gas"].get("scalar", "mass").time)[1:]
        reconstructed = self.d.loader.fine_output_times[1:len(time) + 1]
        rel_diff = (reconstructed - time) / time
        self.assertTrue(all(rel_diff < 1e-10))

    def test_data_dir(self):
        self.assertEqual(self.d.loader.data_dir, code_sample_path)
        #self.assertEqual( self.d.loader.data_dir, code_sample_path + "/out" )

    def test_has_gas(self):
        self.assertTrue("gas" in self.d.fluids)

    def test_has_gasdens(self):
        self.d.fluids["gas"].get("2d", "mass density", 1)

    def test_has_gasprs(self):
        if self.d.loader.eos == 'IDEAL':
            self.d.fluids["gas"].get("2d", "pressure", 1)
        else:
            print(self.eos, 'eos does not support pressure.')

    def test_has_gastemp(self):
        if self.d.loader.eos == 'IDEAL':
            self.d.fluids["gas"].get("2d", "temperature", 1)
        else:
            print(self.eos, 'eos does not support temperature.')

    def test_can_calc_scale_height(self):
        if self.d.loader.eos == 'IDEAL':
            data = self.d.fluids["gas"].get("2d", "pressure scale height", 0)
            #r[36] = 5.194 au
            self.assertTrue(close(data.data[36][0].to('au'),
                                  0.05 * 5.2 * u.au))

        else:
            print(self.eos, 'eos does not support temperature.')

    def test_get_fluid_gas(self):
        self.d.get_fluid("gas")

    def test_gasdens(self):
        #test file structure: density. Also check one specific value
        rho = self.d.fluids["gas"].get("2d", "mass density", 1)
        self.assertEqual(rho.data.shape, (128, 80))
        self.assertTrue(
            close(rho.data[1, 1], 1157507.6839857362 * u.g / u.cm**2))

    def test_gastemp(self):
        #test file structure: temperature. Also check one specific value
        # (checks prs/dens too, if temperature is not explicitly output.)
        if self.d.loader.eos == 'IDEAL':
            tmp = self.d.fluids["gas"].get("2d", "temperature", 0)
            self.assertEqual(tmp.data.shape, (128, 80))
            self.assertTrue(close(tmp.data[36, 0], 121.5 * u.K))
        else:
            print(self.eos, 'eos does not support temperature.')

    def test_gasvrad(self):
        #test file structure: radial velocity
        vrad = self.d.fluids["gas"].get("2d", "vrad", 1)
        self.assertEqual(vrad.data.shape, (128, 80))
        self.assertEqual(len(vrad.grid.get_coordinates("r")), 128)
        self.assertEqual(vrad.data.shape,
                         (len(vrad.grid.get_coordinates("r")),
                          len(vrad.grid.get_coordinates("phi"))))

    def test_output_time(self):
        #test if this output file corresponds to time = 2 Jupiter orbits
        pressure = self.d.fluids["gas"].get("2d", "pressure", 1)
        self.assertTrue(
            close(pressure.time.decompose(), 5 * 11.86 * 3.154e+7 * u.s))

    def test_scalar_mass(self):
        #test if this data element corresponds to time = 2 Jupiter orbits, as well as its value
        mass = self.d.fluids["gas"].get("scalar", "mass")
        self.assertTrue(
            close(mass.data[2].decompose(),
                  2.3811188377819384e+01 * 2e+33 * u.g))
        self.assertTrue(
            close(
                mass.time[2].decompose(),
                2 * 2 * np.pi * np.sqrt((5.2 * 1.49597871e13)**3 * u.cm**3 /
                                        (const.G.cgs * 2e+33 * u.g)).to(u.s)))


def close(x, y):
    if np.abs(x - y) < np.abs(1e-2 * x):
        return True
    else:
        print('failed to match: %.3le to %.3le' % (x.value, y.value))
        return False


if __name__ == '__main__':
    unittest.main()
