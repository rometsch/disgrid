""" Functions to load parameters.
"""
import yaml

def load_params(path):
    param_file = path + "/params.yml"
    with open(param_file, "r") as infile:
        data = yaml.safe_load(infile)
    return data