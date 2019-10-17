from simdata import Data
import numpy as np
import matplotlib.pyplot as plt

d = Data(data_dir='/home/alex/Lib/simdata/samples/PLUTO42/out/')
gas = d.fluids['gas']

i = 1

Sig = d.fluids['gas'].get('2d', 'mass density', i)
time = gas.get('scalar', 'mass').time
print(time)
sig  = Sig.data.T

grid = Sig.grid

r = grid.r_c
phi = grid.phi_c

x, y = np.meshgrid(r, phi)

X, Y = x*np.cos(y), x*np.sin(y)

plt.pcolormesh(X,Y, sig)
plt.show()