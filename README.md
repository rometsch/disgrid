## About disgrid

Simdata is a python library which enables access to simulation data of hydrodynamics planet disk interaction simulations.
Its main features are:

+ load data in a code independend way
  - the same line of code to load data for tool A or B
  - lets you focus on what you want to do with your data
+ automatically detecting the type of simulation code
  - no special book-keeping required to produce comparison plots

## Supported Hydrodynamics Codes

+ [FargoCPT](https://github.com/rometsch/fargocpt) : custom version of the [FARGO code](http://fargo.in2p3.fr/-Legacy-archive-) which is used and maintained at the University of Tübingen.
+ [FARGO3D](http://fargo.in2p3.fr/) : GPU enabled version of the FARGO code supporting 3D calculations
+ [PLUTO 4.2, 4.3](http://plutocode.ph.unito.it/) : magneto-hydrodynamics code based on a Riemann solver

## Getting started

You need python3 to use this package!
Clone this repository and run =python3 setup.py install --user= to install the python package.
Navigate to an output directory of a simulation and enter a python shell and run:

```
from disgrid import Data
# load data from the current directory
d = Data()
# you can also specify a path
d = Data(path_to_data)
# show content of the data
d.avail()
# load some data
mass = d.get(fluid="gas", var="mass", dim="scalar")
# plot the data
plt.plot(mass.time, mass.data)
```

This will print the fluids (gas or dust species) and the planets which disgrid was able to find in your output directory.
