
#ifndef _RIEMANN_H_
#define _RIEMANN_H_

#include "structure.H"

void setRiemannParams( const struct domain * theDomain );
void riemann( struct cell * cL , struct cell * cR, 
              const double r , const double dA , const double dt );


#endif