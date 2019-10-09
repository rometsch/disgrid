import os
import sys
import importlib

class UnknownCodeError(Exception):
    pass

class MultipleCodeError(Exception):
    pass

# functions required for a valid loader module
required_functions = ["identify"]

# store all available loader modules
available = {}

for name in os.listdir(os.path.dirname(os.path.abspath(__file__))):
    if name.endswith(".py"):
        if name in ["__init__.py", "interface.py"]:
            continue
        module_name = name[:-3]
        try:
            module = importlib.import_module(__package__ + "." + module_name)
            missing_functions = [ fname for fname in required_functions if fname not in dir(module) ]
            if len(missing_functions) > 0:
                    print("Warning: module '{}' doesn't supply some require function ({}), ignoring this module".format(module_name, missing_functions), file = sys.stderr)
            else:
                available[module.code_info] = module
        except ImportError:
            print("Warning: module '{}' couldn't be imported".format(module_name), file = sys.stderr)
            raise

def identify_code(path):
    code_list = []
    for key, mod in available.items():
        if mod.identify(path):
            code_list.append(key)
    if len(code_list) == 0:
        raise UnknownCodeError("No known code matches data in '{}'".format(path))
    elif len(code_list) > 1:
        raise MultipleCodeError("Multiple codes identified the data in '{}' which where '{}'".format(path, code_list))
    return code_list[0]

def get_loader(path, loader, **kwargs):
    if loader:
        if "Loader" in dir(loader):
            # assume its a module
            code = loader.code_info
            loader = loader.Loader(path, **kwargs)
        else:
            # assume its a loader object
            loader = loader
            code = type(loader).code_info
    else:
        code = identify_code(path)
        loader = available[code].Loader(path, **kwargs)
    return code, loader
