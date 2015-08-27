
#include "../structure.H"
#include "initial_conditions.H"
#include "ejecta_ICs.H"


void Ejecta_ICs::initial( double * prim , double r )
{

    //The following is consistent with an initial
    //time t = r0/vmax = .05477

    double E = 1.0;
    double M = 1.0;

    double r0 = 0.01;
    double rho0 = 1.0;
    double Pmin = 1e-5;

    double rho = rho0;
    double v   = 0.0;
    double Z   = 0.0;

    double V = 4./3.*M_PI*r0*r0*r0;
    double vmax = sqrt(10./3.*E/M);

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
