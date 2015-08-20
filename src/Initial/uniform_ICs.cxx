
#include "../structure.H"
#include "initial_conditions.H"
#include "uniform_ICs.h"


void Uniform_ICs::initial( double * prim , double r )
{
    prim[RHO] = 1.0;
    prim[PPP] = 1.0;
    prim[VRR] = 0.0;
    prim[ZZZ] = 0.0;
    prim[AAA] = 0.0;
}
