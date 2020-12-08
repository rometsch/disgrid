#!/usr/bin/env python3
import sys, os
import unittest


def get_repo_abspath():
    # add src directory to path for imports
    this_files_dir = os.path.dirname(os.path.abspath(__file__))
    repo_dir = os.path.dirname(this_files_dir)
    return repo_dir


#sys.path.append(os.path.join(get_repo_abspath(), "src"))
sys.path = sys.path[:1] + [os.path.join(get_repo_abspath(), "src")
                           ] + sys.path[1:]

from test_data import *
from test_field import *
from test_fluid import *
from test_grid import *
from test_scalar import *
from test_alias import *
from test_loader_fargocpt import *
from test_loader_fargo3d import *
from test_loader_PLUTO42 import *
from test_loader_PLUTO43 import *

if __name__ == "__main__":
    unittest.main()
