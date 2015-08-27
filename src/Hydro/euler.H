
#ifndef _EULER_H_
#define _EULER_H_

extern "C" {
#include <grackle.h>
}

void setHydroParams( struct domain * theDomain );
void prim2cons( double * prim , double * cons , double dV );
void cons2prim( double * cons , double * prim , double dV , bool verbose=false );
void getUstar( double * prim , double * Ustar , double Sk , double Ss );
void flux( double * prim , double * flux );
void source( double * prim , double * cons , double * grad , 
             double rp , double rm , double dV , double dt , 
             code_units cooling_units, int With_Cooling);
void vel( double * prim1 , double * prim2 , 
          double * Sl , double * Sr , double * Ss );
double mindt( double * prim , double w , double r , double dr );


#endif 