from __future__ import print_function, division

from astropy import units
from astropy import constants as const

hbar     = const.hbar.cgs.value
# # k_b and m_proton need to remain hard-coded
# # in order to stay consistent with the c constants
k_b      = 1.380649e-16 
m_proton = 1.672622e-24 

pc       = units.pc.to(units.cm)
yr       = units.yr.to(units.s)
M_solar  = units.MsolMass.to(units.g)
gamma    = 5/3

metallicity_solar = .02

E_0      = 1e51 # [ergs]

def print_constants():
	"""Prints what the constants SHOULD be,
	in both the python and c files

	This function doesn't actually set anything though.
	That needs to be done by hand.
	"""

	print("G =     ", str("{0:e}").format(const.G.cgs.value))
	print("k_B =   ", str("{0:e}").format(const.k_B.cgs.value))
	print("m_p =   ", str("{0:e}").format(const.m_p.cgs.value))
	print("au =    ", str("{0:e}").format(const.au.cgs.value))
	print("pc =    ", str("{0:e}").format(const.pc.cgs.value))
	print("M_sun = ", str("{0:e}").format(const.M_sun.cgs.value))
	print("R_sun = ", str("{0:e}").format(const.R_sun.cgs.value))
	print("yr =    ", str("{0:e}").format(units.year.to('s')))

	return
