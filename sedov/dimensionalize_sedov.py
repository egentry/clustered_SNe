from __future__ import print_function, division
import numpy as np
import os
from astropy import constants as const 

## Boilerplate path hack to give access to full SNe package
import sys, os
if __package__ is None:
    if os.pardir not in sys.path[0]:
        file_dir = os.path.dirname(__file__)
        sys.path.insert(0, os.path.join(file_dir, os.pardir, os.pardir))

from SNe.constants import hbar, k_b, m_proton, pc, yr, gamma, \
                          metallicity_solar

nstep = 1000
raw_filename = 'spherical_standard_omega0p00_nstep_' + str(nstep).zfill(5) + '.dat'
# make that file name relative to this file
file_dir     = os.path.dirname(__file__)
raw_filename = os.path.join(file_dir, raw_filename)

if not os.path.exists(raw_filename):
    cwd = os.getcwd()
    os.chdir(file_dir)
    os.system("./sedov3 " + str(nstep))
    os.chdir(cwd)

i, x, den, energy, pressure, velocity, cs = np.loadtxt(raw_filename, skiprows=2, unpack=True)

background_density     = m_proton
background_temperature = 1e4


def dimensionalized_sedov(time, metallicity=metallicity_solar, 
                          background_density=background_density,
                          background_temperature=background_temperature):
    """ 
        Scale sedov solution at a given time

        INPUTS:
            - time - [seconds]
    """

    # set mass fractions, to agree with Thornton et al 1998 (i.e. constant Y=.23)
    Z  = metallicity
    Y  = .23
    X  = 1 - Y - Z 
    mu = (2*X + .75*Y + .5*Z)**-1

    rho_0   = background_density
    T_0     = background_temperature
      
    E_0     = 1e51 # erg
    xi_0    = 1.
    r_0     = xi_0 * (E_0 / rho_0)**.2 * (time)**.4
    t_0     = time

        
    # ADD DIMENSIONS

    #  see Eq. 5.23-35 of http://www.astronomy.ohio-state.edu/~ryden/ast825/ch5-6.pdf
    #       except notation will be mixed up -- I can't find a better source at the moment
    r         = r_0 * x
    shock_loc = np.argmax(den)
    r_shock   = r[shock_loc]
    u         = ( (4) / (5 * (gamma+1)) ) * r_shock / t_0                 * velocity / np.max(velocity)
    rho       = ( (gamma+1) / (gamma-1) ) * rho_0                         * den      / np.max(den)
    P         = ( (8) / (25 * (gamma+1))) * rho_0 * (r_shock**2 / t_0**2) * pressure / np.max(pressure)

    v    = (4*np.pi / 3) * (r**3 - np.append([0], r[:-1])**3)
    mass = rho * v

    ener = (1 / (gamma-1)) * P / rho
    c_ad = (gamma * P / rho)**.5
    T    = P / rho * (mu * m_proton / k_b)

    P[(P==0)] = rho[(P==0)] / (mu * m_proton) * k_b * T_0
    T[(T==0)] = T_0

    s = 2.5 - np.log( (rho/(mu * m_proton)) * (2*np.pi*hbar**2 / (mu * m_proton * k_b * T))**1.5 )

    return (r, u, rho, T, c_ad, ener, P, s, mass)


# sedov = dimensionalized_sedov(100 * yr)
