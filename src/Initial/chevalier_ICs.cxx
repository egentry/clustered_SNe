
#include <cmath>
#include <string>

#include "../structure.H"
#include "initial_conditions.H"
#include "chevalier_ICs.H"

const std::string Chevalier_ICs::class_name = "chevalier";

Chevalier_ICs::Chevalier_ICs() : Initial_Conditions( class_name )
{}

int Chevalier_ICs::setICparams( struct domain * theDomain ,
                                const Mass_Loss * mass_loss )
{
    this->set_output_prefix( theDomain );
    this->set_times( theDomain );


    t = theDomain->theParList.t_min;
    return 0;

}

void Chevalier_ICs::initial( double * prim , double r )
{

    const double s = 2.0;
    const double n = 7.0;

    const double g = 1.0;
    const double q = 1.0;

          double rho1 = std::pow(r/t/g,-n)*std::pow(t,-3.);
    const double rho2 = q*std::pow(r,-s);

    const double R0 = std::pow( std::pow(t,n-3.)*std::pow(g,n)/q , 1./(n-s) );
    const double r1 = 0.065*R0;
    if( r<r1 ) rho1 = std::pow(r1/t/g,-n)*std::pow(t,-3.);

    const double rho = rho1+rho2;
    const double v   = (r/t)*rho1/(rho1+rho2);
    const double Z   = rho1/(rho1+rho2);

    const double Pmin = rho*1e-5;

    prim[RHO] = rho;
    prim[PPP] = Pmin;
    prim[VRR] = v;
    prim[ZZZ] = Z;

}
