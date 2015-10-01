
#ifndef _ASCII_H_
#define _ASCII_H_

#include <vector>
#include <string>
#include <cstring> // std::size_t


#include "../structure.H"
#include "../mass_loss.H"
#include "../cooling.H"

void create_checkpoint( struct domain * theDomain , const char * filestart , 
                        const double t );
void overviews(  struct domain * theDomain ,
                Mass_Loss * mass_loss ,
                Cooling * cooling );
void main_overview(  struct domain * theDomain ,
                    Mass_Loss * mass_loss ,
                    Cooling * cooling );
void SNe_overview( struct domain * theDomain );
void inputs_overview( struct domain * theDomain ,
                      Mass_Loss * mass_loss ,
                      Cooling * cooling ); 


std::vector<supernova> read_SNe( const std::string filename);
std::size_t count_lines_in_file( const std::string filename );


#endif