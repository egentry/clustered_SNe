
#include <math.h>

#include "../structure.h"
#include "../constants.h"

double get_dV( double , double );
 
static double mu;     // mean molecular weight -- this is the only time it's used
static double Gamma;  // adiabatic index

static double V_blast;  // volume of the initial blast (1 cell)
static double R_blast;  // outermost boundary of initial blast (1 cell)

static double background_density;      // [g cm^-3]
static double background_temperature;  // [K]



void setICparams( struct domain * theDomain ){
   
   Gamma = theDomain->theParList.Adiabatic_Index;

   const double Y = .23; // helium fraction
   const double Z = theDomain->metallicity; // metals fraction
   const double X = 1 - Y - Z; // hydrogen mass fraction

   mu = 1. / (2*X + .75*Y + .5*Z); // mean molecular weight

   int i;
   int n_blast = 2;
   V_blast = 0.;
   R_blast = 0.;
   for( i=0 ; i<n_blast ; ++i )
   {
      struct cell * c = &(theDomain->theCells[i]);
      double rp = c->riph;
      double rm = rp - c->dr; 
      V_blast += get_dV( rp , rm );
      R_blast = rp;
   }

   background_density     = theDomain->background_density;
   background_temperature = theDomain->background_temperature;

}

void initial( double * prim , double r ){
   double rho,Pp;

   const double M_sun = 1.989100e+33;  // [g]

   if( r <= R_blast ){
      // inner guard cell will also be set to these values
      //   but that's okay, since guard cell needs to have same density + pressure
      // this shouldn't affect energy conservation
      const double M_blast = 3 * M_sun;       // [g]
      const double E_blast = 1e51 / M_blast;  // [erg g^-1]

      rho = M_blast / V_blast;
      Pp  = rho * (Gamma-1) * E_blast;
   }else{
      // initial background conditions
      rho = background_density;
      const double T_0 = background_temperature;  // [K]
      Pp  = (rho / (mu * m_proton)) * k_boltzmann * T_0;
   }
   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;
}
