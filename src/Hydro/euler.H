
#ifndef _EULER_H_
#define _EULER_H_

extern "C" {
#include <grackle.h>
}

#include "../cooling.H"

void setHydroParams( const struct domain * theDomain );
void prim2cons( const double * prim , double * cons , const double dV );
void cons2prim( const double * cons , double * prim , 
                const double dV , const bool verbose=false );
double E_int_from_cons( const struct cell * c );
double E_int_from_prim( const struct cell * c , 
                        const double gamma );
void getUstar( const double * prim , double * Ustar , 
               const double Sk , const double Ss );
void flux( const double * prim , double * flux );
void source( const double * prim , double * cons , const double * grad , 
                   double * dE_cool , 
             const double rp , const double rm , 
             const double dV , const double dt , 
             Cooling * cooling );
void conduction( struct cell * cL , struct cell * cR, 
              const double dA , const double dt );
void vel( const double * prim1 , const double * prim2 , 
          double * Sl , double * Sr , double * Ss );
double mindt( const double * prim , const double w , 
              const double r , const double dr );


#endif 