# a loader object orchestrates the loading of data from files
# it 
# it provides a framework for caching data and provides function which return field objects
class Interface:
    def __init__(self, path):
        self.path = path
        self.fluids = {}
        self.vectors = {}
        self.particlegroups = {}
        self.planets = {}
        self.parameters = {}
        self.meta = {}

    def scout(self):
        # find all variables
        # adjust self.fluids
        # adjust self.field_loaders
        pass

    def get(self, *args, **kwargs):
        pass

    def get_output_time(self, n):
        pass
