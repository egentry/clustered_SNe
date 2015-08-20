
#include "../structure.H"
#include "initial_conditions.H"
#include "isentropic_ICs.h"

void Isentropic_ICs::initial( double * prim , double r )
{
    double R2 = r*r;
    prim[RHO] = 1.0 + 3.0*exp(-80.*R2);
    prim[PPP] = pow(prim[RHO],5./3.);
    prim[VRR] = 0.0;
    prim[ZZZ] = 0.0;
    prim[AAA] = 0.0;
}