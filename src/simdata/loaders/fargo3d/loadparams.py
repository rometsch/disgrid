""" Functions to load parameters from fargo3d output files using the summmary files.
"""
import os
import re

def find_first_summary(dataDir):
    return "summary{}.dat".format(find_first_summary_number(dataDir))


def find_first_summary_number(dataDir):
    return find_summary_numbers(dataDir)[0]


def find_summary_numbers(dataDir):
    ptrn = re.compile(r"summary(\d+).dat")
    summaries = []
    for f in os.listdir(dataDir):
        m = re.search(ptrn, f)
        if m:
            n = int(m.groups()[0])
            summaries.append(n)
    summaries.sort()
    return summaries


def getParamFromSummary(dataDir, param):
    return getParamsFromNthSummary(
        dataDir, find_first_summary_number(dataDir))[param.lower()]


def getParamsFromNthSummary(dataDir, n):
    # parse the Nth summary file to get all
    search_active = False
    parameters = {}
    with open(os.path.join(dataDir, "summary{}.dat".format(n))) as f:
        for line in f:
            line = line.strip()
            if not search_active:
                # look for the parameter section identifier
                if line == "PARAMETERS SECTION:":
                    search_active = True
                continue
            if line == "" or line[0] in ["#", "="]:
                continue
            if line.startswith("*** Input file: "):
                parameters["config path"] = line.split(":")[-1].strip()
                break
            parts = [s.strip() for s in line.split()]
            try:
                val = int(parts[1])
            except ValueError:
                try:
                    val = float(parts[1])
                except ValueError:
                    val = parts[1]
            parameters[parts[0].lower()] = val
    return parameters
