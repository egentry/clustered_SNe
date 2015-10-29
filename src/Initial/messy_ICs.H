
#ifndef _MESSY_ICS_H_
#define _MESSY_ICS_H_

#include <string>

#include "initial_conditions.H"
#include "../structure.H"

class Messy_ICs : public Initial_Conditions
{

public:

    Messy_ICs();

	void initial( double * prim , double r );

    static const std::string class_name;

    void possibly_extend_grid( struct domain * theDomain);
    double find_shock( const struct domain * theDomain ) const;

private:

};

#endif
