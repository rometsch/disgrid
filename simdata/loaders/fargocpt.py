code_info = ( "fargocpt", "0.1", "testloader")

import os
from . import interface

def identify(path):
    return "fargo" in os.listdir(path)

class Loader(interface.LoaderInterface):
    
    def __init__(self, path):
        self.grid = None

class FieldQuantity2d(interface.FieldQuantity):

    def load_data(self, n):
        pass
        
class FieldRho(interface.FieldQuantity):

    def load_data(self, n):
        pass
    
    
