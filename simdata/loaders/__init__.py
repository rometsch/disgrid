import os
import sys
import importlib

# functions required for a valid loader module
required_functions = ["identify"]

# store all available loader modules
available = {}


for name in os.listdir(os.path.dirname(os.path.abspath(__file__))):
    if name.endswith(".py"):
        if name == "__init__.py":
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
        

