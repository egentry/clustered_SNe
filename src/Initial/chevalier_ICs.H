
#ifndef _CHEVALIER_ICS_H_
#define _CHEVALIER_ICS_H_

#include <string>

#include "initial_conditions.H"
#include "../structure.H"
#include "../mass_loss.H"

class Chevalier_ICs : public Initial_Conditions
{

public:

    Chevalier_ICs();

    int setICparams( struct domain * theDomain ,
                     const Mass_Loss * mass_loss );
    void initial( double * prim , double r );

    static const std::string class_name;

    void possibly_extend_grid( struct domain * theDomain);
    double find_shock( const struct domain * theDomain ) const;

private:

    double t;

};

#endif