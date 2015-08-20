
#ifndef _EJECTA_ICS_H_
#define _EJECTA_ICS_H_

#include "initial_conditions.H"
#include "../structure.H"

class Ejecta_ICs : public Initial_Conditions
{

public:

    void initial( double * prim , double r );

private:

};

#endif