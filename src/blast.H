
#ifndef _blast_H_
#define _blast_H_

#include <vector>

#include "structure.H"

int add_single_blast( struct domain * theDomain , const double E_blast = 1e51 );
int add_blasts( struct domain * theDomain );

int get_SNe( const double , std::vector<double>&SNe_times ,
                      const double ,
                      const unsigned int );

#endif