""" Functions to load units.
"""
import yaml
import astropy

def load_units(path):
    param_file = path + "/units.yml"
    with open(param_file, "r") as infile:
        data = yaml.safe_load(infile)

    units = {}
    units["mass"] = astropy.units.Unit(data["M"])
    units["time"] = astropy.units.Unit(data["T"])
    units["length"] = astropy.units.Unit(data["L"])

    return units


def get_unit_from_powers(unitpowers, units):
    unit = 1.0
    for k, p in unitpowers.items():
        unit = unit * units[k]**p
    return unit