
#ifndef _EULER_H_
#define _EULER_H_

extern "C" {
#include <grackle.h>
}

void setHydroParams( const struct domain * theDomain );
void prim2cons( const double * prim , double * cons , const double dV );
void cons2prim( const double * cons , double * prim , 
                const double dV , const bool verbose=false );
void getUstar( const double * prim , double * Ustar , 
               const double Sk , const double Ss );
void flux( const double * prim , double * flux );
void source( const double * prim , double * cons , const double * grad , 
             const double rp , const double rm , 
             const double dV , const double dt , 
             code_units cooling_units, const int With_Cooling);
void vel( const double * prim1 , const double * prim2 , 
          double * Sl , double * Sr , double * Ss );
double mindt( const double * prim , const double w , 
              const double r , const double dr );


#endif 