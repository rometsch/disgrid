from .data import Data

try:
    from .smurf_plugin import SmurfData as SData
except (ModuleNotFoundError,ImportError):
    pass

try:
    from .server_plugin import NData
except ImportError:
    pass