class LoaderInterface:

    def __init__(self):
        pass

class Quantity:

    def __init__(self, loader):
        self.loader = loader
        self.dataDir = loader.dataDir

class FieldQuantity(Quantity):
    pass
