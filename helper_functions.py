from __future__ import division

import numpy as np

from numba import jit, njit

from constants import m_proton, hbar, k_b, gamma


@jit
def calculate_mean_molecular_weight(metallicity):
    Z = metallicity
    Y = .23
    X = 1 - Y - Z
    mu = (2*X + .75*Y + .5*Z)**-1
    return mu 

@jit
def calculate_mass(density, dV):
    mass = density * dV
    return mass

@jit
def calculate_kinetic_energy(mass, velocity):
    E_kin = (1/2) * mass * velocity**2
    return E_kin.sum()

@jit
def calculate_internal_energy(mass, pressure, density):
    E_int = mass * (1/(gamma-1)) * pressure / density
    return E_int.sum()

@jit
def calculate_momentum(mass, velocity):
    momentum = mass * velocity
    return momentum.sum()

@jit
def calculate_c_ad(pressure, density):
    c_ad = np.sqrt(gamma * pressure / density)
    return c_ad

@jit
def calculate_entropy(temperature, density, mu):
    s = 2.5 - np.log( (density/(mu * m_proton)) \
                      * (2*np.pi*hbar**2 / (mu * m_proton * k_b * temperature))**1.5 )
    return s

@jit
def calculate_temperature(pressure, density, mu):
    temperature = pressure / density * (mu * m_proton / k_b)
    return temperature