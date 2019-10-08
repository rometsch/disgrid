#!/usr/bin/env python3
import os
from setuptools import setup

setup(  name="simdata"
        ,version="0.1"
        ,description="scripts to hydro simulation output"
        ,author="Thomas Rometsch, Peter Rodenkirch"
        ,packages=["simdata", "simdata.loaders"]
        ,install_requires=["numpy>=1.17", "astropy"]
        )
