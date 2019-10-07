import unittest
from simdata import data
from simdata.loaders import fargo3dmultifluid
import astropy.units as u
import astropy.constants as const
import numpy as np

code_sample_path = "samples/fargo3dmultifluid"

class TestFargo3DMultifluidLoader(unittest.TestCase):

    def setUp(self):
        self.d = data.Data(code_sample_path)

    def test_identify_code_via_data(self):
        self.assertEqual( self.d.loader.code_info, ( "fargo3d", "2.0", "multifluid") )

    def test_identify_code_directly(self):
        self.assertTrue( fargo3dmultifluid.identify(code_sample_path) )

    def test_units(self):
        self.assertEqual( self.d.loader.units['length'], 5.2 * 1.49597871e13 * u.cm)
        self.assertEqual( self.d.loader.units['mass'] , 1.9891e+33 * u.g)
        self.assertEqual( self.d.loader.units['time'],  np.sqrt((5.2 * 1.49597871e13)**3 * u.cm**3 / (const.G.cgs * 1.9891e+33 * u.g)).to(u.s))

    def test_data_dir(self):
        self.assertEqual( self.d.loader.data_dir, code_sample_path + "/outputs" )

    def test_has_gas(self):
        self.assertTrue("gas" in self.d.fluids)

    def test_has_gasdens(self):
        self.d.fluids["gas"].get("2d", "mass density", 0)

    def test_get_fluid_gas(self):
        self.d.get_fluid("gas")

    def test_get_fluid_dust(self):
        self.d.get_fluid("dust1")
        self.d.get_fluid("dust2")
        self.d.get_fluid("dust3")

    def test_gasdens(self):
        rho = self.d.fluids["gas"].get("2d", "mass density", 2)
        self.assertEqual( rho.data.shape , (64, 32) )
        self.assertEqual( rho.data[1,1] , 4.285694082764542 * u.g / u.cm**2 )

    def test_dustdens(self):
        rho = self.d.fluids["dust1"].get("2d", "mass density", 2)
        self.assertEqual( rho.data.shape , (64, 32) )
        self.assertAlmostEqual( rho.data[1,1].value , 9.190435006078763e-07 )

    def test_gasvrad(self):
        vrad = self.d.fluids["gas"].get("2d", "vrad", 2)
        self.assertEqual( vrad.data.shape , (64, 32) )
        self.assertEqual( len(vrad.grid.get_coordinates("r")) , 64 )
        self.assertEqual( vrad.data.shape, (len(vrad.grid.get_coordinates("r")), len(vrad.grid.get_coordinates("phi")) ) )

    def test_output_time(self):
        energy = self.d.fluids["gas"].get("2d", "energy density", 2)
        self.assertTrue( np.abs(energy.time.decompose() - 12.5663706200000007*5.9551995752415031e+07*u.s) < 1e-3*energy.time.decompose())

    # def test_vector_mass(self):
    #     mass = self.d.fluids["gas"].get("vector", "mass")
    #     self.assertEqual(mass.data[2].decompose(), 7.8098717547161137e-03*1.9889199999999999e+33*u.g)
    #     self.assertEqual(mass.time[2].decompose(), 1.2566370620000000e-01*5.9551995752415031e+07*u.s)

    # def test_multidim_vector(self):
    #     ekin = self.d.fluids["gas"].get("vector", "kinetic energy")
    #     self.assertEqual( ekin.get(5).shape, (2,))
    #     self.assertEqual( ekin.get(5)[0],  ekin.get(5, axis="r") )
    #     self.assertEqual( ekin.get(5)[1],  ekin.get(5, axis="phi") )

    # def test_vector_time(self):
    #     ekin = self.d.fluids["gas"].get("vector", "kinetic energy")
    #     t, v = ekin.get(slice(2,4), return_time=True)
    #     self.assertTrue( all(t.decompose() == [1.2566370620000000e-01*5.9551995752415031e+07, 1.8849555930000000e-01*5.9551995752415031e+07]*u.s) )

    # def test_particlegroups(self):
    #     self.assertTrue( "planets" in self.d.particlegroups )
    #     self.assertTrue( "omega frame" in self.d.particlegroups["planets"].variable_loaders )

    # def test_single_planets_via_particlegroups(self):
    #     x = self.d.particlegroups["planets"].get("x", 0)
    #     #self.assertEqual( x , None)
    #     self.assertEqual( x[1].decompose(), 0.600009762159255167*7.7790892764000000e+13*u.cm)

    # def test_multiple_planets_via_particlegroups(self):
    #     x = self.d.particlegroups["planets"].get("x")
    #     self.assertEqual( x.data.shape, (2,201))
    #     #self.assertEqual( x , None)
    #     self.assertTrue( all (x[1].decompose() == np.array([0.600009762159255167, 0.997363963701585199])*7.7790892764000000e+13*u.cm ))

    # def test_multiple_planets_multi_axes_via_particlegroups(self):
    #     pos = self.d.particlegroups["planets"].get("position")
    #     self.assertEqual( pos.data.shape, (2,201,2))
    #     self.assertTrue( all (pos[:,1,0].decompose() == np.array([0.600009762159255167, 0.997363963701585199])*7.7790892764000000e+13*u.cm ))
    #     vel = self.d.particlegroups["planets"].get("velocity")
    #     self.assertEqual( vel.data.shape, (2,201,2))
    #     self.assertTrue( all (vel[:,1,0].decompose() == np.array([0.000306735295912537458, 0.0718294173695262217])*1.3062684429152035e+06*u.cm/u.s ))

    def test_planet(self):
        x = self.d.planets[0].get("x")
        #self.assertEqual( x , None)
        self.assertTrue( close(x[1].decompose(), 9.22999995481561974* 5.2 * 1.49597871e13 * u.cm ))

    # def test_planet_multi_axes(self):
    #     pos = self.d.planets[0].get("position")
    #     self.assertEqual( pos.data.shape, (201,2) )
    #     self.assertEqual( pos[1,0].decompose(), 0.600009762159255167*7.7790892764000000e+13*u.cm )


def close(x,y):
    return np.abs(x-y)<np.abs(1e-5*x)
    
if __name__ == '__main__':
    unittest.main()
