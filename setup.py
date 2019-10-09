#!/usr/bin/env python3
import os
from setuptools import setup, find_namespace_packages

setup(  name="simdata"
        ,version="0.1"
        ,description="scripts to load hydro simulation output"
        ,author="Thomas Rometsch, Peter Rodenkirch"
        ,package_dir={'': 'src'}
        ,packages=find_namespace_packages(where="src")
        ,install_requires=["numpy>=1.17", "astropy"]
        )
