import unittest
from simdata import data
from simdata.loaders import PLUTO42
import astropy.units as u
import astropy.constants as const
import numpy as np
import os

from run import get_repo_abspath
code_sample_path = os.path.join(get_repo_abspath(),"samples/PLUTO42/out")

class TestPLUTO42Loader(unittest.TestCase):

    def setUp(self):
        self.d = data.Data(code_sample_path)

    def test_identify_code_via_data(self):
        self.assertEqual( self.d.loader.code_info, ( "PLUTO", "4.2", "vanilla") )

    def test_identify_code_directly(self):
        self.assertTrue( PLUTO42.identify(code_sample_path) )

    def test_units(self):
        self.assertTrue( close( 1.*self.d.loader.units['length'], 5.2 * 1.49597871e13 * u.cm) )
        self.assertTrue( close( 1.*self.d.loader.units['mass'], 2e33 * u.g) )
        self.assertTrue( close( 1.*self.d.loader.units['time'], 2*np.pi*np.sqrt((5.2 * 1.49597871e13)**3 * u.cm**3 / (const.G.cgs * 2e33 * u.g)).to(u.s)) )

    #def test_parameters(self):
    #    self.assertEqual( self.d.parameters["amin"] , 0.001)
    #    self.assertEqual( self.d.parameters["planetconfig"] , 'setups/hd163296/hd163296.cfg')

    def test_fine_output_times(self):
        #initial time is 0 -> removed first element
        time = (self.d.fluids["gas"].get("scalar", "mass").time)[1:]
        reconstructed = self.d.loader.fine_output_times[1:len(time)+1]
        rel_diff = (reconstructed - time)/time
        self.assertTrue( all( rel_diff < 1e-10 ) )

    def test_data_dir(self):
        self.assertEqual( self.d.loader.data_dir, code_sample_path )
        #self.assertEqual( self.d.loader.data_dir, code_sample_path + "/out" )

    def test_has_gas(self):
        self.assertTrue("gas" in self.d.fluids)

    def test_has_gasdens(self):
        self.d.fluids["gas"].get("2d", "mass density", 1)

    def test_get_fluid_gas(self):
        self.d.get_fluid("gas")

    #def test_get_fluid_dust(self):
        #self.d.get_fluid("dust1")
        #self.d.get_fluid("dust2")
        #self.d.get_fluid("dust3")

    def test_gasdens(self):
        #test density file structure: density. Also check one specific value
        rho = self.d.fluids["gas"].get("2d", "mass density", 1)
        self.assertEqual( rho.data.shape , (128, 80) )
        self.assertTrue( close( rho.data[1,1], 1157507.6839857362 * u.g / u.cm**2) )

    #def test_dustdens(self):
    #    rho = self.d.fluids["dust1"].get("2d", "mass density", 2)
    #    self.assertEqual( rho.data.shape , (64, 32) )
    #    self.assertAlmostEqual( rho.data[1,1].value , 9.190435006078763e-07 )

    def test_gasvrad(self):
        #test file structure: radial velocity
        vrad = self.d.fluids["gas"].get("2d", "vrad", 1)
        self.assertEqual( vrad.data.shape , (128, 80) )
        self.assertEqual( len(vrad.grid.get_coordinates("r")) , 128 )
        self.assertEqual( vrad.data.shape, (len(vrad.grid.get_coordinates("r")), len(vrad.grid.get_coordinates("phi")) ) )

    def test_output_time(self):
        #test if this output file corresponds to time = 2 Jupiter orbits
        pressure = self.d.fluids["gas"].get("2d", "pressure", 1)
        self.assertTrue( close( pressure.time.decompose(), 5*11.86*3.154e+7*u.s) )

    def test_scalar_mass(self):
        #test if this data element corresponds to time = 2 Jupiter orbits, as well as its value
        mass = self.d.fluids["gas"].get("scalar", "mass")
        self.assertTrue( close(mass.data[2].decompose(), 2.3811188377819384e+01*2e+33 * u.g) )
        self.assertTrue( close(mass.time[2].decompose(), 2*2*np.pi*np.sqrt((5.2 * 1.49597871e13)**3 * u.cm**3 / (const.G.cgs * 2e+33 * u.g)).to(u.s)) )

    #def test_scalar_torque_planet_1(self):
    #    torq = self.d.fluids["dust1"].get("scalar", "torque planet 1")
    #    T = np.sqrt((5.2 * 1.49597871e13)**3 * u.cm**3 / (const.G.cgs * 1.9891e+33 * u.g)).to(u.s)
    #    L = 5.2 * 1.49597871e13 * u.cm
    #    M = 1.9891e+33 * u.g
    #    self.assertTrue(close(torq.data[3].decompose().cgs, -8.46861771884e-11*M*L**2/T**2))
    #    self.assertTrue(close(torq.time[3].decompose().cgs, 2.51327412288*T))

    #def test_1d_torque_planet_1(self):
    #    torq = self.d.fluids["dust2"].get("1d", "torque planet 1", 1)
    #    T = np.sqrt((5.2 * 1.49597871e13)**3 * u.cm**3 / (const.G.cgs * 1.9891e+33 * u.g)).to(u.s)
    #    L = 5.2 * 1.49597871e13 * u.cm
    #    M = 1.9891e+33 * u.g
    #    self.assertTrue(close(torq.data[15].decompose().cgs, -8.30347786497707e-19*M*L**2/T**2))


    # def test_multidim_scalar(self):
    #     ekin = self.d.fluids["gas"].get("scalar", "kinetic energy")
    #     self.assertEqual( ekin.get(5).shape, (2,))
    #     self.assertEqual( ekin.get(5)[0],  ekin.get(5, axis="r") )
    #     self.assertEqual( ekin.get(5)[1],  ekin.get(5, axis="phi") )

    # def test_scalar_time(self):
    #     ekin = self.d.fluids["gas"].get("scalary", "kinetic energy")
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

    #def test_planet(self):
    #    x = self.d.planets[0].get("x")
    #    #self.assertEqual( x , None)
    #    self.assertTrue( close(x[1].decompose(), 9.22999995481561974* 5.2 * 1.49597871e13 * u.cm ))

    # def test_planet_multi_axes(self):
    #     pos = self.d.planets[0].get("position")
    #     self.assertEqual( pos.data.shape, (201,2) )
    #     self.assertEqual( pos[1,0].decompose(), 0.600009762159255167*7.7790892764000000e+13*u.cm )


def close(x,y):
    return np.abs(x-y)<np.abs(1e-2*x)

if __name__ == '__main__':
    unittest.main()
