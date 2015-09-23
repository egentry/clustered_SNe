
#ifndef _MISC_H_
#define _MISC_H_

extern "C" {
#include <grackle.h>
}

#include "cooling.H"
#include "structure.H"

double getmindt( const struct domain * theDomain );
double get_mean_molecular_weight( const double Z );
void set_wcell( struct domain * theDomain );
void adjust_RK_cons( struct domain * theDomain , const double RK );
void move_cells( struct domain * theDomain , const double dt);
void calc_dr( struct domain * theDomain );
void fix_negative_energies( struct domain * theDomain );
void calc_prim( struct domain * theDomain );
void radial_flux( struct domain * theDomain , const double dt );
void add_source( struct domain * theDomain , const double dt ,
                 Cooling * cooling );
void longandshort( const struct domain * theDomain , 
                   double * L , double * S , 
                   int * iL , int * iS );
void AMR( struct domain * theDomain );


#endif
