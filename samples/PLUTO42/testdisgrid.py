from disgrid import Data
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

d = Data('/home/alex/Lib/disgrid/samples/PLUTO42/out/')
#d = Data('/home/alex/Lib/disgrid/samples/PLUTO43/out/')
print(d.code)
gas = d.fluids['gas']

i = 1

Dens0 = d.fluids['gas'].get('2d', 'mass density', 0)
Dens = d.fluids['gas'].get('2d', 'mass density', i)

print(Dens.data.shape)
print(d.loader.geometry)
time = gas.get('scalar', 'mass').time
#print(time)
dens  = Dens.data.T
dens0 = Dens0.data.T

grid = Dens.grid

r = grid.r_c
phi = grid.phi_c

if '4.3' in d.code:
	
	theta = grid.theta_c
	fig, ax = plt.subplots(1,2, figsize=(8,3))
	ax[0].set_aspect(1)

	x, y = np.meshgrid(r, phi)
	X, Y = x*np.cos(y), x*np.sin(y)
	midplane = (dens/dens0)[:,0,:].decompose().cgs.value

	img = ax[0].pcolormesh(X,Y, np.log10(midplane))
	cb = plt.colorbar(img, ax=ax[0])
	cb.set_label(r'$\log \rho/\rho_0$')

	x, y = np.meshgrid(r, np.pi/2*u.radian-theta)
	X, Y = x*np.cos(y), x*np.sin(y)
	meridian = (dens/dens0).mean(axis=0).decompose().cgs.value

	img = ax[1].pcolormesh(X,Y, np.log10(meridian))
	img = ax[1].pcolormesh(X,-Y, np.log10(meridian))
	cb = plt.colorbar(img, ax=ax[1])
	cb.set_label(r'$\log \rho/\rho_0$')
	
	for axis in ax: axis.set_xlabel('x [5.2 au]')
	ax[0].set_ylabel('y [5.2 au]')
	ax[1].set_ylabel('z [5.2 au]')

else:
	
	fig, ax = plt.subplots(1,1, figsize=(5,4))
	ax.set_aspect(1)

	x, y = np.meshgrid(r, phi)
	X, Y = x*np.cos(y), x*np.sin(y)
	midplane = (dens/dens0).decompose().cgs.value

	img = ax.pcolormesh(X,Y, np.log10(midplane))
	cb = plt.colorbar(img, ax=ax)
	cb.set_label(r'$\log \Sigma/\Sigma_0$')

	ax.set_xlabel('x [5.2 au]')
	ax.set_ylabel('y [5.2 au]')
	
plt.tight_layout()
plt.show()