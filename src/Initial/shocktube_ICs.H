
#ifndef _SHOCKTUBE_ICS_H_
#define _SHOCKTUBE_ICS_H_

#include <string>

#include "initial_conditions.H"
#include "../structure.H"

class Shocktube_ICs : public Initial_Conditions
{

public:

    Shocktube_ICs();

    void initial( double * prim , double r );

    static const std::string class_name;

private:

};

#endif