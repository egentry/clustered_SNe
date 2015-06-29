#ifndef _CONSTANTS_H
#define _CONSTANTS_H


#define gamma  (5./3)      // adiabatic index
#define mu  .67      // mean molecular weight -- should probably be calculated properly using metallicity


// PHYSICAL CONSTANTS SET USING ASTROPY VALUES
// but when using python, you should also use these values, due to round-off error
#define G  6.673840e-08
#define k_boltzmann 1.380649e-16
#define m_proton  1.672622e-24
#define comp  (k_boltzmann / (m_proton * mu))

#define c_V  ((1./(gamma - 1)) * k_boltzmann / (mu * m_proton)) // specific heat capacity, const. volume

#define AU  1.495979e+13
#define pc  3.085678e+18
#define M_sun  1.989100e+33
#define R_sun  6.955080e+10

#define yr 3.155760e+07

// EVERYTHING BELOW HERE SHOULD PROBABLY BE MOVED TO ANOTHER FILE


// set these more properly later
#define M_SN         (3 * M_sun)
#define E_SN         (1e51 / M_SN) // Energy per unit mass

#define R_total      (1e2 * pc ) // total size, initially
#define V_0          (1. / (m_proton)) // initial background density
#define T_0          1e4             // initial background temperature

#define V_core       (1. / (100 * m_proton))
#define R_core       0 // radius of innermost zone edge (unresolved) 
#define M_core       ((4./3) * M_PI * pow(R_core,3) * V_core) // mass internal to R_core
#define P_core       0
#define Q_core       0
#define E_core       0
#define T_core       0
#define C_ad_core    sqrt(gamma * P_core * V_core)

// #define V_background      = V_0
// #define T_background      = T_0

#define R_inner      (1. * pc)   // size of beginning of sedov phase
#define V_inner      ((4./3) * M_PI * (pow(R_inner,3) - pow(R_core,3)) / M_SN)
#define T_inner      (E_SN / c_V)
#define P_inner      (comp * T_inner / V_inner)

#endif	