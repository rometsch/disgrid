import os
from astropy import units as u
from subprocess import run, PIPE

log_file_options = ['pluto.log', '../log.log', 'pluto.0.log']


def loadUnits(dataDir, dimensions, filecache=None):
    # extract lines for units from pluto log file (specify at top of file!)
    unit_lines = []
    in_unit_block = False

    use_log_file = None
    for log_file in log_file_options:
        if os.path.exists(os.path.join(dataDir, log_file)):
            use_log_file = log_file
            break

    if not use_log_file:
        raise NotImplementedError(
            'Cannot find PLUTO log file. If compiled in PLUTO 4.3 serial mode, it will be implemented differently.'
        )

    log_file_path = os.path.join(dataDir, use_log_file)
    if filecache is not None:
        log_file_path = filecache(log_file_path)
    grep_res = run([
        "grep", "-C", "7", "Normalization\ Units",
        log_file_path
    ],
                   stdout=PIPE)
    lines = grep_res.stdout.decode("utf-8").splitlines()[-8:]

    for line in lines:
        line = line.strip()
        if line == "> Normalization Units:":
            in_unit_block = True
            continue
        if in_unit_block:
            if line == "":
                continue
            elif line[0] == ">":
                break
            else:
                unit_lines.append(line)
    # parse lines into units dict
    units = {}
    for line in unit_lines:
        parts = line.split()
        name = parts[0].strip(":[]").lower()
        value = float(parts[1])
        unit_str = parts[2].strip(",()")
        if unit_str == "X":
            unit_str = parts[-1].strip("()")
        unit = u.Unit(unit_str.replace('gr', 'g').replace('sec', 's'))
        units[name] = u.Unit(value * unit)
    # calculate mass unit from density and length
    units["density"] *= units["length"]**(3 - dimensions
                                          )  # fix for non-3d density
    units["mass"] = u.Unit(
        (units["density"] * units["length"]**dimensions).decompose()).cgs
    units["time"] = units["time"] / (1 * units["length"]).value  #why?
    return units
