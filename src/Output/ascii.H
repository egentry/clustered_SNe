
#ifndef _ASCII_H_
#define _ASCII_H_

#include <cstring>

#include "../structure.H"

void output( struct domain * theDomain , const char * filestart , 
             const double t );
int overview( struct domain * theDomain );
std::size_t count_lines_in_file( const std::string filename );


#endif