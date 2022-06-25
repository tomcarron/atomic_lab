"""Script for functions used in the analysis for the FTIR experiment in the atomic and molecular lab course"""
import numpy as np
import pandas as pd
#import math
#import decimal
#from sigfig import round as rnd

""" reduced_mass: Function to calculate the reduced mass of two nuclei, with errors."""


def reduced_mass(m1, dm1, m2, dm2):
    mu = (m1 * m2) / (m1 + m2)
    top = (dm1**2) * (m1**4) + (dm2**2) * (m2**4)
    bottom = (m1 + m2) ** 2
    dmu = np.sqrt(top) / bottom
    return mu, dmu


""" inertia: Calculates the moment of inertia from the equilibrium bond length
Takes Be and error in cm, returns moment of inertia I in amu angstrom^2 and error"""


def inertia(Be, dBe):
    const = 16.8576314
    I = const / Be
    dI = const * I * dBe / Be
    return I, dI


""" B_e: Calculates the equilibrium bond length from the constants B0 and B1
Also calculate the error in the equilibrium bond length from their errors"""


def B_e(B0, dB0, B1, dB1):
    Be = (3 * B0 - B1) * 0.5
    dBe = Be * np.sqrt((1.5 * dB0) ** 2 + (-0.5 * dB1) ** 2)
    return Be, dBe


""" bond_length: Calculates the experimentally derived bond length and error from the
moment of inertia and reduced mass and their errors"""


def bond_length(I, dI, mu, dmu):
    r = np.sqrt(I / mu)
    dr = r * np.sqrt(
        (dI / (2 * np.sqrt(I * mu))) ** 2
        + ((dmu * np.sqrt(I)) / (2 * mu ** (1.5))) ** 2
    )
    return r, dr
