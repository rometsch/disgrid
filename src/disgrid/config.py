import os
import json

home_path = os.path.join(os.path.expanduser("~"), ".disgrid")


class Config:
    def __init__(self):
        self.config_file = os.path.join(home_path, "config.json")
        if not os.path.exists(home_path):
            self.create()
        else:
            self.load()

    def __getitem__(self, key):
        return self.data[key]

    def create(self):
        os.makedirs(home_path)
        self.load()
        self.save()
        
        
    def __setitem__(self, key, val):
        self.data[key] = val
        self.save()

    def save(self):
        self.data["type"] = "disgrid config"
        self.data["version"] = "0.1"
        with open(self.config_file, "w") as out_file:
            out_file.write(json.dumps(self.data, indent=4))

    def load(self):
        try:
            with open(self.config_file, "r") as in_file:
                self.data = json.load(in_file)
        except FileNotFoundError:
            self.data = {}
