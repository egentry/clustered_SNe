#ifndef _COOLING_H
#define _COOLING_H

void   calc_cooling(double V[], double E[], double U[], double dE[], int zones, double metallicity, double dt, code_units my_units);  
code_units   setup_cooling();  


#endif