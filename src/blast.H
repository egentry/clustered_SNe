
#ifndef _blast_H_
#define _blast_H_

#include <vector>
#include <string>

#include "structure.H"
#include "mass_loss.H"


int add_single_blast( struct domain * theDomain , const double M_blast,
                                                  const double M_blast_Z,
                                                  const double E_blast = 1e51);

int add_blasts( struct domain * theDomain );

std::vector<supernova> get_SNe( const double ,
                                const double ,
                                const unsigned int ,
                                const Mass_Loss * mass_loss );
std::vector<supernova> read_SNe(  const std::string filename );


bool sort_by_lifetime( const supernova &a, const supernova &b);






#endif