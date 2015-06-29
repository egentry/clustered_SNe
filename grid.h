#ifndef _GRID_H
#define _GRID_H


	void swap(double **a, double**b);
	void expand_grid(int *zones_ptr, 
	    double **U_old_ptr, double **U_new_ptr,
	    double **R_old_ptr, double **R_new_ptr,
	    double **V_old_ptr, double **V_new_ptr,
	    double **T_old_ptr, double **T_new_ptr,
	    double **E_old_ptr, double **E_new_ptr, double **dE_ptr,
	    double **P_old_ptr, double **P_new_ptr, 
	    double **Q_old_ptr, double **Q_new_ptr,
	    double **H_old_ptr, double **H_new_ptr,
	    double **C_ad_ptr,
	    double **Mass_ptr,  double **M_int_ptr);


	void create_Mass_array(double Mass[], double R[], double V[], int zones);
	void create_M_int_array(double M_int[], double Mass[], int zones);

	void enforce_boundary_conditions(int zones, 
	    double U_old[], double U_new[], 
	    double R_old[], double R_new[], 
	    double V_old[], double V_new[], 
	    double T_old[], double T_new[], 
	    double E_old[], double E_new[], 
	    double P_old[], double P_new[], 
	    double Q_old[], double Q_new[],
	    double H_old[], double H_new[] );

	void enforce_boundary_conditions_initial(int zones, 
	    double U_old[], double U_new[], 
	    double R_old[], double R_new[], 
	    double V_old[], double V_new[], 
	    double T_old[], double T_new[], 
	    double E_old[], double E_new[], 
	    double P_old[], double P_new[], 
	    double Q_old[], double Q_new[],
	    double H_old[], double H_new[] );

	void enforce_boundary_conditions_Dirichlet(int zones, 
	    double U_old[], double U_new[], 
	    double R_old[], double R_new[], 
	    double V_old[], double V_new[], 
	    double T_old[], double T_new[], 
	    double E_old[], double E_new[],
	    double H_old[], double H_new[] );

	void enforce_boundary_conditions_Neumann(int zones, 
	    double P[], 
	    double Q[] );



	double interp(double f_inner, double f_outer);

	double harmonic_mean(double a, double b);


#endif