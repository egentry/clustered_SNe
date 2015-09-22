
#ifndef _EJECTA_ICS_H_
#define _EJECTA_ICS_H_

#include <string>

#include "initial_conditions.H"
#include "../structure.H"

class Ejecta_ICs : public Initial_Conditions
{

public:

    Ejecta_ICs();

    void initial( double * prim , double r );

    static const std::string class_name;
    
private:


};

#endif