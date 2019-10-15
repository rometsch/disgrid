from .data import Data

# import simscripts integration when simscripts is installed
try:
    import simscripts.search
    from .simscriptsdata import SData
except ImportError as e:
    print(e)
    pass
