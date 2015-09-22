
#ifndef _UNIFORM_ICS_H_
#define _UNIFORM_ICS_H_

#include <string>

#include "initial_conditions.H"
#include "../structure.H"

class Uniform_ICs : public Initial_Conditions
{

public:

    Uniform_ICs();

    void initial( double * prim , double r );

    static const std::string class_name;

private:

};

#endif