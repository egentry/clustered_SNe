#ifndef _FLUID_ADIABATIC_H
#define _FLUID_ADIABATIC_H

	double calc_E_kin( double Mass[], double U[], int zones );
	double calc_E_int( double Mass[], double T[], int zones );
	double calc_E_grav(double Mass[], double M_int[], double R[], int zones );


	double calc_Momentum(double Mass[], double U[], int zones);



#endif