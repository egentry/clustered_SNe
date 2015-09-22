
#define _USE_MATH_DEFINES // for M_PI
#include <cmath>
#include <string>

#include "../structure.H"
#include "initial_conditions.H"
#include "ejecta_ICs.H"

const std::string Ejecta_ICs::class_name = "ejecta";

Ejecta_ICs::Ejecta_ICs() : Initial_Conditions( class_name )
{}

void Ejecta_ICs::initial( double * prim , double r )
{

    //The following is consistent with an initial
    //time t = r0/vmax = .05477

    const double E = 1.0;
    const double M = 1.0;

    const double r0 = 0.01;
    const double rho0 = 1.0;
    const double Pmin = 1e-5;

    double rho = rho0;
    double v   = 0.0;
    double Z   = 0.0;

    const double V = 4./3.*M_PI*r0*r0*r0;
    const double vmax = std::sqrt(10./3.*E/M);

    if( r < r0 ){
      rho = M/V;
      v = vmax*r/r0;
      Z = 1.0;
    }

    prim[RHO] = rho;
    prim[PPP] = Pmin;
    prim[VRR] = v;
    prim[ZZZ] = Z;

}
