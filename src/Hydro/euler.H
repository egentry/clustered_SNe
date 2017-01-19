
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
double E_int_from_cons( const double * cons );
double E_int_from_prim( const struct cell * c , 
                        const double gamma );
void getUstar( const double * prim , double * Ustar , 
               const double Sk , const double Ss );
void flux( const double * prim , double * flux );
void source( struct cell * c , 
             const double rp , const double rm , 
             const double dV , const double dt , 
             const double R_shock , Cooling * cooling );
void conduction( struct cell * cL , struct cell * cR, 
              const double dA , const double dt );
void vel( const double * prim1 , const double * prim2 , 
          double * Sl , double * Sr , double * Ss );
double mindt( const double * prim , const double w , 
              const double r , const double dr );

double calc_T(const double * prim);

void calc_multiphase_prim(const struct cell * c, 
    double * prim_hot, double * prim_cold,
    double * V_hot,    double * V_cold );

void subgrid_conduction( struct cell * c , const double dt );


#endif 