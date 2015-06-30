#include <math.h>

#include "../paul.h"
#include "../constants.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r , int i , int Nr , double V_blast){
   double rho,Pp;
   if( i < 2 ){
      rho = M_SN / V_blast;
      Pp  = rho * (gamma-1) * E_SN;
   }else{
      rho = 1.  / V_0;
      Pp  = (rho / (mu * m_proton)) * k_boltzmann * T_0;
   }
   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;
}
