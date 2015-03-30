from __future__ import print_function
import numpy as np

i, x, den, energy, pressure, velocity, cs = np.loadtxt('sedov/spherical_standard_omega0p00.dat', skiprows=2, unpack=True)

gamma   = 5./3 			# adiabatic index

k_b     = 1.38046e-16 	# boltzmann's constant
mu      = .67			# mean molecular weight [dimensionless]
m_p     = 1.6598e-24 	# proton mass [g]
pc      = 3.086e18 		# 1 parsec [cm]
yr 		= 3.154e7    	# year in [s]


rho_0   = m_p
  
r_0     = 1. * pc
E_0     = 1e51 # erg
xi_0    = x[np.argmax(den)]
t_0     = ((r_0 / xi_0)**5 * rho_0 / E_0 )**.5
T_0     = 1e4 # background temperature [K]

print(r_0)

print( "current time: ", t_0, "[s]" )
print( "current time: ", t_0 / yr, "[yr]" ) # print time in years

xi_max  = x[-1]
d_xi    = x[-1] - x[-2]


new_zones = 250 - 1 # total zones MUST BE ODD
x_new           = np.linspace(xi_max + d_xi, xi_max +  d_xi * new_zones, new_zones)
den_new         = np.ones(new_zones) * den[-1]
energy_new      = np.zeros(new_zones)
pressure_new    = np.zeros(new_zones)
velocity_new    = np.zeros(new_zones)
cs_new          = np.zeros(new_zones)

x           = np.append(x,          x_new)
den         = np.append(den,        den_new)
energy      = np.append(energy,     energy_new)
pressure    = np.append(pressure,   pressure_new)
velocity    = np.append(velocity,   velocity_new)
cs          = np.append(cs,         cs_new)

    
# ADD DIMENSIONS

rho = ( (gamma + 1) / (gamma-1) ) * rho_0 * den / np.max(den)
v   = 1. / rho
u   = ( (4) / (5 * (gamma+1)) ) * r_0 / t_0 * velocity / np.max(velocity)
P   = ( (8) / (25 * (gamma+1))) * rho_0 * r_0**2 / t_0**2 * pressure / np.max(pressure)
r   = r_0 * x / xi_0

outside     = np.where(P == 0 )
P[outside]  = rho[outside] / (mu * m_p) * k_b * T_0

ener = (1 / (gamma-1)) * P * v
c_ad = (gamma * P / rho)**.5
T    = P / rho * (mu * m_p / k_b)
q    = np.zeros(T.shape)


np.savetxt("input", np.array([ u, r, v, T, ener, P, q, c_ad]).transpose(), comments='',
    header="# t_0 [s] \n" + str(t_0) + "\n" + "# \t u [cm s^-1] \t\t\t\t r [cm] \t\t\t\t V [cm^3 g^-1] \t\t\t\t T [K] \t\t\t\t E [erg g^-1] \t\t\tP [dyne cm^-2] \t\t\t\t q [dyne cm^-2] \t\t c_ad [cm s^-1]")


