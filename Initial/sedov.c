
#include <math.h>

#include "../structure.h"
#include "../constants.h"

static double Gamma = 5./3;  // adiabatic index

void setICparams( struct domain * theDomain ){
   
   Gamma = theDomain->theParList.Adiabatic_Index;

}

void initial( double * prim , double r , int i , int Nr , double V_blast){
   double rho,Pp;
   if( i < 2 ){
      double M_blast = 3 * M_sun;       // [g]
      double E_blast = 1e51 / M_blast;  // [erg g^-1]

      rho = M_blast / V_blast;
      Pp  = rho * (Gamma-1) * E_blast;
   }else{
      // initial background conditions
      rho = 1. * m_proton;
      double T_0 = 1e4;  // [K]
      Pp  = (rho / (mu * m_proton)) * k_boltzmann * T_0;
   }
   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;
}
