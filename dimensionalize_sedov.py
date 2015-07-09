from __future__ import print_function
import numpy as np
import os
from astropy import constants as const 

nstep = 1000
raw_filename = 'sedov/spherical_standard_omega0p00_nstep_' + str(nstep).zfill(4) + '.dat'

if not os.path.exists(raw_filename):
    os.chdir("sedov")
    os.system("./sedov3 " + str(nstep))
    os.chdir("..")

i, x, den, energy, pressure, velocity, cs = np.loadtxt(raw_filename, skiprows=2, unpack=True)

# h_bar does not need to be in agreement with constants.h
# but it should agree with the value being used in visualize.ipynb
hbar = const.hbar.cgs.value

k_b     = 1.380649e-16  # boltzmann's constant
m_p     = 1.672622e-24  # proton mass [g]
pc      = 3.085678e+18  # 1 parsec [cm]
yr      = 3.155760e+07  # year in [s]

metallicity_solar      = .02 #
background_density     = m_p
background_temperature = 1e4


def main(time, metallicity=metallicity_solar, 
               background_density=background_density,
               background_temperature=background_temperature):
    """ 
        Scale sedov solution at a given time

        INPUTS:
            - time - [seconds]
    """

    # set mass fractions, to agree with Thornton et al 1998 (i.e. constant Y=.23)
    Z = metallicity
    Y = .23
    X = 1 - Y - Z 
    mu = (2*X + .75*Y + .5*Z)**-1.

    # THESE UNITS MUST AGREE WITH THE CONSTANTS.H FILE
    # Alternative they can be generated using get_constants.py
    gamma   = 5./3          # adiabatic index



    rho_0   = background_density
    T_0     = background_temperature
      
    E_0     = 1e51 # erg
    xi_0    = 1.
    r_0     = xi_0 * (E_0 / rho_0)**.2 * (time)**.4
    t_0     = time


        
    # ADD DIMENSIONS

    #  see Eq. 5.23-35 of http://www.astronomy.ohio-state.edu/~ryden/ast825/ch5-6.pdf
    #       except notation will be mixed up -- I can't find a better source at the moment
    r       = r_0 * x
    r_shock = r[np.argmax(den)]
    u       = ( (4.) / (5 * (gamma+1)) )  * r_shock / t_0                 * velocity / np.max(velocity)
    rho     = ( (gamma+1) / (gamma-1) )   * rho_0                         * den      / np.max(den)
    P       = ( (8.) / (25. * (gamma+1))) * rho_0 * (r_shock**2 / t_0**2) * pressure / np.max(pressure)

    v    = (4*np.pi / 3) * (r**3 - np.append([0,0], r[:-2])**3)
    mass = rho * v

    ener = (1 / (gamma-1)) * P / rho
    c_ad = (gamma * P / rho)**.5
    T    = P / rho * (mu * m_p / k_b)

    P[(P==0)] = rho[(P==0)] / (mu * m_p) * k_b * T_0
    T[(T==0)] = T_0

    # print(x[np.argmax(pressure)])
    s = 2.5 - np.log( (rho/(mu * m_p)) * (2*np.pi*hbar**2 / (mu * m_p * k_b * T))**1.5 )

    return (r, u, rho, T, c_ad, ener, P, s, mass)


yr = 3.155760e+07
# main(16.4614780406 * yr)