
#ifndef _ASCII_H_
#define _ASCII_H_

#include <vector>
#include <string>
#include <cstring> // std::size_t


#include "../structure.H"
#include "../mass_loss.H"
#include "../cooling.H"

void output( struct domain * theDomain , const char * filestart , 
             const double t );
int overview( struct domain * theDomain ,
              Mass_Loss * mass_loss ,
              Cooling * cooling );
std::vector<supernova> read_SNe( const std::string filename);
std::size_t count_lines_in_file( const std::string filename );


#endif