
#ifndef _CHEVALIER_ICS_H_
#define _CHEVALIER_ICS_H_

#include "initial_conditions.H"
#include "../structure.H"
#include "../mass_loss.H"

class Chevalier_ICs : public Initial_Conditions
{

public:

    int setICparams( struct domain * theDomain ,
                     const Mass_Loss * mass_loss );
    void initial( double * prim , double r );

private:

    double t;

};

#endif