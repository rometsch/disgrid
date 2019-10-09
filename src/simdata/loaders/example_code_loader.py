code_info = ( "exmaple_code", "0.1", "test")

import os
from . import interface

def identify(path):
    try:
        with open(os.path.join(path, "text.txt")) as f:
            for line in f:
                if line.strip() == "very specific content!":
                    return True
            return False
    except FileNotFoundError:
        return False

class Loader(interface.Interface):
    pass
