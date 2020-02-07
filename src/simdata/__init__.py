from .data import Data

try:
    from .smurf_plugin import RemoteData
    from .smurf_plugin import SmurfData as SData
except ImportError:
    pass
