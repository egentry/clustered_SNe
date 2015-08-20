
#ifndef _COOLING_H_
#define _COOLING_H_

extern "C" {
#include <grackle.h>
}

double calc_cooling( double * prim ,  double * cons , 
                     double dt , code_units my_units );

code_units setup_cooling( struct domain * theDomain );                

#endif 
