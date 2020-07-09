""" Functions to grid data from fargo3d output files.
"""
import os

import numpy as np

from . import loadparams


def loadRadius(dataDir, unit, interface=False):
    r = np.genfromtxt(os.path.join(dataDir, 'domain_y.dat')) * unit
    r = r[3:-3]  # remove ghost cells
    dr = r[1:] - r[:-1]
    if not interface:
        r = 0.5 * (r[1:] + r[:-1])
    return (r, dr)


def loadPhi(dataDir, interface=False):
    #phiMin, phiMax, Nphi = np.genfromtxt(os.path.join(dataDir, 'dimensions.dat'), usecols=(0,1,6))
    #phi = np.linspace(phiMin, phiMax, Nphi)
    phi = np.genfromtxt(os.path.join(dataDir, 'domain_x.dat'))
    if not interface:
        phi = 0.5 * (phi[1:] + phi[:-1])
    return phi


def loadMeshGrid(dataDir, unit):
    # return a meshgrid for the disk to plot data
    R, Phi = loadMeshGridPolar(dataDir, unit)
    X = R * np.cos(Phi)
    Y = R * np.sin(Phi)
    return (X, Y)


def loadMeshGridPolar(dataDir, unit):
    phi = loadPhi(dataDir)
    r, dr = loadRadius(dataDir, unit)
    Phi, R = np.meshgrid(phi, r)
    return (R, Phi)


def loadNcells(dataDir):
    # Nphi, Nr = np.genfromtxt(os.path.join(dataDir, 'dimensions.dat'), usecols=(6,7), dtype=int)
    Nphi = int(loadparams.getParamFromSummary(dataDir, "Nx"))
    Nr = int(loadparams.getParamFromSummary(dataDir, "Ny"))
    return (Nphi, Nr)
