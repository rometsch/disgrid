class UnknownCodeError(Exception):
    pass

class MultipleCodeError(Exception):
    pass

from . import example_loader
from . import fargo3d
from . import fargo3d_1_3
from . import fargocpt
from . import fargocpt_1_0
from . import fargocpt_1_1
from . import fargocpt_1_2
from . import PLUTO42
from . import PLUTO43

available = {}
available[example_loader.code_info] = example_loader
available[fargo3d.code_info] = fargo3d
available[fargo3d_1_3.code_info] = fargo3d_1_3
available[fargocpt.code_info] = fargocpt
available[fargocpt_1_0.code_info] = fargocpt_1_0
available[fargocpt_1_1.code_info] = fargocpt_1_1
available[fargocpt_1_2.code_info] = fargocpt_1_2
available[PLUTO42.code_info] = PLUTO42
available[PLUTO43.code_info] = PLUTO43


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
            if key[0] == loader[0] and key[1] == loader[1]:
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
