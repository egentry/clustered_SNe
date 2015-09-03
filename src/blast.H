
#ifndef _blast_H_
#define _blast_H_

#include <vector>
#include <string>

#include "structure.H"



void add_winds( struct domain * theDomain, const double dt );
void add_winds_constant( struct domain * theDomain , const double dt, 
                  const double wind_velocity = 1e3 * 1000 * 100 );

int add_single_blast( struct domain * theDomain , const double M_blast,
                                                  const double M_blast_Z,
                                                  const double E_blast = 1e51);

int add_blasts( struct domain * theDomain );

std::vector<supernova> get_SNe( const double ,
                                const double ,
                                const unsigned int );
std::vector<supernova> read_SNe(  const std::string filename );


bool sort_by_lifetime( const supernova &a, const supernova &b);
double get_ejecta_mass( const double M_initial );
double get_ejecta_mass_Z( const double M_initial );





#endif