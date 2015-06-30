#ifndef _ICS_H
#define _ICS_H


	double sedov(int zones, double U[], double R[], double V[], double T[], 
	            double E[], double P[], double Q[], double H[], double C_ad[]);

	double new_blast(int zones,
	    double U[], double R[], double V[], double T[],
	    double E[], double P[], double Q[], double H[], double C_ad[]);

	double new_blast_spread(int zones,
	    double U[], double R[], double V[], double T[],
	    double E[], double P[], double Q[], double H[], double C_ad[]);	

	double restart(int *zones_ptr, int *k_ptr,
		double U[], double R[], double V[], double T[],
	    double E[], double P[], double Q[], double H[], double C_ad[],
	    char mod_filename[]);

#endif