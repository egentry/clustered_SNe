from __future__ import print_function
import numpy as np
import os
import sys

gamma   = 5./3 			# adiabatic index

k_b     = 1.38046e-16 	# boltzmann's constant
mu      = .67			# mean molecular weight [dimensionless]
m_p     = 1.6598e-24 	# proton mass [g]
pc      = 3.086e18 		# 1 parsec [cm]
yr 		= 3.154e7    	# year in [s]

rho_0   = m_p
  
r_0     = 1. * pc
r_outer = 100 * pc
E_0     = 1e50 # erg


sedov_zones = 160
if len(sys.argv) > 1:
	sedov_zones = int(sys.argv[1])

outer_zones = sedov_zones*4
if (sedov_zones + outer_zones)%2 is 0:
	outer_zones += 1   # ensure odd number of total zones

if len(sys.argv) > 2:
	outer_zones = int(sys.argv[2])
	if (sedov_zones + outer_zones)%2 is 0:
		raise ValueError("sedov_zones + outer_zones must give odd number of total zones")

outer_by_radius = False
if len(sys.argv) > 3:
	r_outer = float(sys.argv[3])
	outer_by_radius = True



raw_filename = "sedov/spherical_standard_omega0p00_nstep_" + str(sedov_zones).zfill(4) + ".dat"

if not os.path.isfile(raw_filename):
	os.chdir("sedov")
	os.system("./sedov3 " + str(sedov_zones))
	os.chdir("..")

i, x, den, energy, pressure, velocity, cs = np.loadtxt(raw_filename, skiprows=2, unpack=True)


xi_0    = x[np.argmax(den)]
t_0     = ((r_0 / xi_0)**5 * rho_0 / E_0 )**.5
T_0     = 1e4 # background temperature [K]

print("r_0:", r_0)

print( "current time: ", t_0, "[s]" )
print( "current time: ", t_0 / yr, "[yr]" ) # print time in years

xi_max  = x[-1]
if outer_by_radius is True:
	d_xi = (xi_0 / r_0) * (r_outer - r_0) / outer_zones
else:
	d_xi = x[-1] - x[-2]

if (sedov_zones + outer_zones)%2 is 0:
	outer_zones += 1   # total zones MUST BE ODD
x_new       = np.linspace(xi_max + d_xi, xi_max +  d_xi * outer_zones, outer_zones)

den_new         = np.ones(outer_zones) * den[-1]
energy_new      = np.zeros(outer_zones)
pressure_new    = np.zeros(outer_zones)
velocity_new    = np.zeros(outer_zones)
cs_new          = np.zeros(outer_zones)

x           = np.append(x,          x_new)
den         = np.append(den,        den_new)
energy      = np.append(energy,     energy_new)
pressure    = np.append(pressure,   pressure_new)
velocity    = np.append(velocity,   velocity_new)
cs          = np.append(cs,         cs_new)

# ADD DIMENSIONS

r   = r_0 * x / xi_0
rho = ( (gamma + 1) / (gamma-1) ) * rho_0 * den / np.max(den)
v   = 1. / rho
u   = ( (4) / (5 * (gamma+1)) ) * r_0 / t_0 * velocity / np.max(velocity)
P   = ( (8) / (25 * (gamma+1))) * rho_0 * r_0**2 / t_0**2 * pressure / np.max(pressure)

outside     = np.where(P == 0 )
P[outside]  = rho[outside] / (mu * m_p) * k_b * T_0

ener = (1 / (gamma-1)) * P * v
c_ad = (gamma * P / rho)**.5
T    = P / rho * (mu * m_p / k_b)
q    = np.zeros(T.shape)


np.savetxt("input", np.array([ u, r, v, T, ener, P, q, c_ad]).transpose(), comments='',
    header="# t_0 [s] \n" + str(t_0) + "\n" + "# \t u [cm s^-1] \t\t\t\t r [cm] \t\t\t\t V [cm^3 g^-1] \t\t\t\t T [K] \t\t\t\t E [erg g^-1] \t\t\tP [dyne cm^-2] \t\t\t\t q [dyne cm^-2] \t\t c_ad [cm s^-1]")


