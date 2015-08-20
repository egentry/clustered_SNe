
#ifndef _TIMESTEP_H_
#define _TIMESTEP_H_

#include "structure.H"

void substep( struct domain * theDomain , double RK , 
              double dt , int first_step , int last_step );
void timestep( struct domain * theDomain , double dt );


#endif