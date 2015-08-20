
#ifndef _CHEVALIER_ICS_H_
#define _CHEVALIER_ICS_H_

#include "initial_conditions.H"
#include "../structure.H"

class Chevalier_ICs : public Initial_Conditions
{

public:

    int setICparams( struct domain * theDomain );
    void initial( double * prim , double r );

private:

    double t;

};

#endif