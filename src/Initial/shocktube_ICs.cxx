
#include <string>

#include "../structure.H"
#include "initial_conditions.H"
#include "shocktube_ICs.H"

const std::string Shocktube_ICs::class_name = "shocktube";

Shocktube_ICs::Shocktube_ICs() : Initial_Conditions( class_name )
{}

void Shocktube_ICs::initial( double * prim , double r )
{
    double rho,Pp;
    if( r < 0.25 )
    {
        rho = 1.0;
        Pp  = 1.0;
    }
    else
    {
        rho = 0.1;
        Pp  = 0.1;
    }
    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[VRR] = 0.0;
    prim[ZZZ] = 0.0;
}
