
#ifndef _SHOCKTUBE_ICS_H_
#define _SHOCKTUBE_ICS_H_

#include "initial_conditions.H"
#include "../structure.H"

class Shocktube_ICs : public Initial_Conditions
{

public:

    void initial( double * prim , double r );

private:

};

#endif