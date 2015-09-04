
#ifndef _COOLING_H_
#define _COOLING_H_

extern "C" {
#include <grackle.h>
}

double calc_cooling( const double * prim ,  const double * cons , 
                     const double dt , code_units my_units );

code_units setup_cooling( const struct domain * theDomain );                

#endif 
