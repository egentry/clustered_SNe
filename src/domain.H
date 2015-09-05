
#ifndef _domain_H_
#define _domain_H_

#include "structure.H"
#include "Initial/initial_conditions.H"

int  setupDomain( struct domain * theDomain , 
                  Initial_Conditions * ICs ,
                  Mass_Loss * mass_loss );
void freeDomain( struct domain * theDomain );
void check_dt( struct domain * theDomain , double * dt );
void possiblyOutput( struct domain * theDomain , int override );


#endif