
#include <cmath>

#include "../structure.H"
#include "initial_conditions.H"
#include "messy_ICs.H"


void Messy_ICs::initial( double * prim , double r )
{

    const double r0 = 1.0;//1.3941e16;
    const double t0 = 1.0;//86400*45.;
    const double rho_csm = 1.0;//1.3177e-17;

    const double rc = 0.1*r0;
    const double rt = 0.31654*rc;

    const double day = t0/45.;
    const double t_ej = 5.*day;
    const double rho_ej  = 3.0307*rho_csm;

    const double v0 = r0/t0;
    const double Poverrho = 1e-5*v0*v0;

    double v   = 0.0;
    double Z   = 0.0;
    double rho = rho_csm;

    if( r < rc ){
        rho = rho_ej*std::pow(rc/r,10.);
        v = r/t_ej;
        Z = 1.0;
        if( r < rt ) rho = rho_ej*std::pow(rc/rt,10.)*rt/r;
    }
    if( r>2.*rc ){
        rho *= .0001;
    }

    prim[RHO] = rho;
    prim[PPP] = rho*Poverrho;
    prim[VRR] = v;
    prim[ZZZ] = Z;

}
