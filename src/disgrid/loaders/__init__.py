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

_loader_src_dir = os.path.dirname(os.path.abspath(__file__))

def import_loader(module_str):
    try:
        module = importlib.import_module(module_str)
        missing_functions = [
            fname for fname in required_functions
            if fname not in dir(module)
        ]
        if len(missing_functions) > 0:
            print(
                "Warning: module '{}' doesn't supply some require function ({}), ignoring this module"
                .format(module_name, missing_functions),
                file=sys.stderr)
            return None
        else:
            return module
    except ImportError:
        print("Warning: module '{}' couldn't be imported".format(
            module_name),
              file=sys.stderr)
        raise

for name in os.listdir(_loader_src_dir):
    if name.endswith(".py"):
        if name in ["__init__.py", "interface.py"]:
            continue
        module_name = name[:-3]
        module_str = __package__ + "." + module_name
        mod = import_loader(module_str)
        if mod is not None:
            available[mod.code_info] = mod
    abs_path = os.path.join(_loader_src_dir, name)
    if os.path.isdir(abs_path):
        if "loader.py" in os.listdir(abs_path):
            module_str = __package__ + "." + name
            mod = import_loader(module_str)
            if mod is not None:
                available[mod.code_info] = mod

def identify_code(path, choices=available):
    code_list = []
    for key, mod in choices.items():
        if mod.identify(path):
            code_list.append(key)
    if len(code_list) == 0:
        raise UnknownCodeError(
            "No known code matches data in '{}'".format(path))
    elif len(code_list) > 1:
        raise MultipleCodeError(
            "Multiple codes identified the data in '{}' which where '{}'".
            format(path, code_list))
    return code_list[0]


def get_loader(path, loader, **kwargs):
    """
    Get a loader for the simulation data inside path.

    Parameters
    ----------
    path : str
        Path to the data.
    loader : str or :obj:`disgrid.loaders.interface.Interface`
        None or hint (str) or acutal loader for the data.

    Returns
    -------
    disgrid.loaders.interface.Interface
        Loader object to access data.
    """
    code = None
    choices = available
    # check for direct specification of code
    if isinstance(loader, (tuple, list)):
        for key, mod in available.items():
            if key == tuple(loader):
                code = key
                try:
                    loader = mod.Loader(path, **kwargs)
                except Exception:
                    code = None
                break
    # check for simulation code names hints
    elif isinstance(loader, str):
        choices = {
            key: available[key]
            for key in available if loader in key[0]
        }
        if len(choices) == 0:
            choices = available
    # handle loader objects/classes
    elif loader is not None:
        if "Loader" in dir(loader):
            # assume its a module
            code = loader.code_info
            loader = loader.Loader(path, **kwargs)
        else:
            # assume its a loader object
            loader = loader
            code = type(loader).code_info
    # if no hints or loader was given, test all available
    if code is None:
        code = identify_code(path, choices)
        loader = available[code].Loader(path, **kwargs)
    return code, loader
