from collections import OrderedDict
import importlib
import sys

class UnknownCodeError(Exception):
    pass

class MultipleCodeError(Exception):
    pass

# library of the available loader implementations
loader_lib = {
    # ( "name", "version", "description" ) : module
    ("fargo3d", "1.3", "custom") : "fargo3d_v1_3_custom"
    }

required_functions = ["identify"]

# attempt to import all available loader modules and store them
# into the loader_lib dict, remove the loader entry if the
# file fails to import
to_remove = []
for key, modulestr in loader_lib.items():
    try:
        loader_lib[key] = importlib.import_module(module_str)
    except ImportError:
        print("Warning: module '{}' could not be imported, its deleted from the loader list".format(module_str), file = sys.stderr)
        to_remove.append(key)
for k in to_remove:
    del loader_lib[k]

to_remove = []
# check the loader implementations for required functions
for key, mod in loader_lib.items():
    for fname in require_functions:
    if not fname in dir(mod):
        print("Warning: module '{}' does not supply require function, its deleted from the loader list".format(module_str), file = sys.stderr)
        if not key in to_remove:
            to_remove.append(key)
for k in to_remove:
    del loader_lib[k]
del to_remove

def identify_code(path):
    code_list = []
    for key, mod in loader_lib.values():
        if mod.identify_code(path):
            code_list.append(code)
    if len(code_list) == 0:
        raise UnknownCodeError("No known code matches data in '{}'".format(path))
    elif len(code_list) > 1:
        raise MultipleCodeError("Multiple codes identified the data in '{}' which where '{}'".format(path, code_list))
    return code_list[0]

def get_loader(code):
    if code == "fargo3d_v1.3":
        import fargo3d_v1_3.loader as loader

class Simdata:
    def __init__(self, path, **kwargs):
        self.path = path
        self.code = self.identify_code()
        self.loader = loader_lib[code].Loader(path, **kwargs)
        self.fields = {}
        self.nbodysystems = {}
        self.parameters = {}
        self.meta = {}

    
