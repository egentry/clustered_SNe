
#ifndef _ISENTROPIC_ICS_H_
#define _ISENTROPIC_ICS_H_

#include "initial_conditions.H"
#include "../structure.H"

class Isentropic_ICs : public Initial_Conditions
{

public:

    void initial( double * prim , double r );

private:

};

#endif