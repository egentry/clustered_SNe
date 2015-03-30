from __future__ import print_function
import numpy as np

def main(time):
	""" 
		Scale sedov solution at a given time

		INPUTS:
			- time - [seconds]
	"""

	i, x, den, energy, pressure, velocity, cs = np.loadtxt('sedov/spherical_standard_omega0p00.dat', skiprows=2, unpack=True)

	gamma   = 5./3 			# adiabatic index

	k_b     = 1.38046e-16 	# boltzmann's constant
	mu      = .67			# mean molecular weight [dimensionless]
	m_p     = 1.6598e-24 	# proton mass [g]
	pc      = 3.086e18 		# 1 parsec [cm]
	yr 		= 3.154e7    	# year in [s]

	rho_0   = m_p
	  
	E_0     = 1e51 # erg
	xi_0    = 1.115
	r_0     = xi_0 * (E_0 / rho_0)**.2 * (time)**.4
	t_0     = time

	T_0     = 1e4 # background temperature [K]

	xi_max  = x[-1]
	d_xi    = x[-1] - x[-2]

	    
	# ADD DIMENSIONS

	r   = r_0 * x / xi_0
	u   = ( (4) / (5 * (gamma+1)) ) * r_0 / t_0 * velocity / np.max(velocity)
	rho = ( (gamma + 1) / (gamma-1) ) * rho_0 * den / np.max(den)
	P   = ( (8) / (25 * (gamma+1))) * rho_0 * r_0**2 / t_0**2 * pressure / np.max(pressure)

	v 	 = (4*np.pi / 3) * (r**3 - np.append([0,0], r[:-2])**3)
	# print(r[0])
	# print(v)
	mass = rho * v


	ener = (1 / (gamma-1)) * P / rho
	c_ad = (gamma * P / rho)**.5
	T    = P / rho * (mu * m_p / k_b)

	P[(P==0)] = rho[(P==0)] / (mu * m_p) * k_b * T_0

	s = np.log(P * np.power(rho, -gamma))

	print(r_0)

	return (r, u, rho, T, c_ad, ener, P, s, mass)

	# np.savetxt("input", np.array([ u, r, v, T, ener, P, q, c_ad]).transpose(), comments='',
	#     header="# t_0 [s] \n" + str(t_0) + "\n" + "# \t u [cm s^-1] \t\t\t\t r [cm] \t\t\t\t V [cm^3 g^-1] \t\t\t\t T [K] \t\t\t\t E [erg g^-1] \t\t\tP [dyne cm^-2] \t\t\t\t q [dyne cm^-2] \t\t c_ad [cm s^-1]")


yr = 3.154e7
main(16.4614780406 * yr)