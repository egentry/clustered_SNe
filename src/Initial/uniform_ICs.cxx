
#include <string>

#include "../structure.H"
#include "initial_conditions.H"
#include "uniform_ICs.H"


const std::string Uniform_ICs::class_name = "uniform";

Uniform_ICs::Uniform_ICs() : Initial_Conditions( class_name )
{}

void Uniform_ICs::initial( double * prim , double r )
{
    prim[RHO] = 1.0;
    prim[PPP] = 1.0;
    prim[VRR] = 0.0;
    prim[ZZZ] = 0.0;
}
