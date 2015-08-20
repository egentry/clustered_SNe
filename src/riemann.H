
#ifndef _RIEMANN_H_
#define _RIEMANN_H_

#include "structure.H"

void setRiemannParams( struct domain * theDomain );
void riemann( struct cell * cL , struct cell * cR, 
              double r , double dA , double dt );


#endif