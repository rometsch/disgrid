#!/usr/bin/env python3
import os
from distutils.core import setup

setup(  name="simdata"
        ,version="alpha 1"
        ,description="scripts to hydro simulation output"
        ,author="Thomas Rometsch, Peter Rodenkirch"
        ,packages=["simdata", "simdata.loaders"]
        )
