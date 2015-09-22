
#ifndef _TIMESTEP_H_
#define _TIMESTEP_H_

#include "cooling.H"
#include "structure.H"
#include "Initial/initial_conditions.H"


void substep( struct domain * theDomain , double RK , 
              double dt , int first_step , int last_step,
              Initial_Conditions * ICs ,
              Cooling * cooling );
void timestep( struct domain * theDomain , const double dt,
              Initial_Conditions * ICs ,
              Cooling * cooling );


#endif