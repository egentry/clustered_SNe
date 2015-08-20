
#ifndef _MISC_H_
#define _MISC_H_

extern "C" {
#include <grackle.h>
}

#include "structure.H"

double getmindt( struct domain * theDomain );
void set_wcell( struct domain * theDomain );
void adjust_RK_cons( struct domain * theDomain , double RK );
void move_cells( struct domain * theDomain , double RK , double dt);
void calc_dr( struct domain * theDomain );
void fix_negative_energies( struct domain * theDomain );
void calc_prim( struct domain * theDomain );
void radial_flux( struct domain * theDomain , double dt );
void add_source( struct domain * theDomain , double dt );
void longandshort( struct domain * theDomain , 
                   double * L , double * S , 
                   int * iL , int * iS );
void AMR( struct domain * theDomain );


#endif
