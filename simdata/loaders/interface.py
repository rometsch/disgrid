# a loader object orchestrates the loading of data from files
# it 
# it provides a framework for caching data and provides function which return field objects
class Loader:
    def __init__(self, path):
        self.path = path
        self.

class Quantity:
    def __init__(self, loader):
        self.loader = loader
        self.dataDir = loader.dataDir

class Field3d(Quantity):
    pass

class Field2d(Quantity):
    pass
        
class Field1d(Quantity):
    pass

class Scalar(Quantity):
    pass
