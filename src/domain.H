
#ifndef _domain_H_
#define _domain_H_

int setupDomain( struct domain * theDomain );
void setupCells( struct domain * theDomain );
void freeDomain( struct domain * theDomain );
void check_dt( struct domain * theDomain , double * dt );
void possiblyOutput( struct domain * theDomain , int override );


#endif