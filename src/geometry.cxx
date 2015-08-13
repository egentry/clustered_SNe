
#include <math.h>

#include "structure.H"

double get_dA( double r ){
   return( 4.*M_PI*r*r );
}

double get_moment_arm( double rp , double rm ){

   // ============================================= //
   //
   //  A relatively arbitrary way to get a cell centered radius
   //
   //  Inputs:
   //    - rp  - outer cell boundary
   //    - rm  - inner cell boundary
   //
   //  Returns:
   //    - r   - cell "center"
   //
   //  Side effects:
   //    None
   //
   //  Notes:
   //    - I don't see any physical motivation behind this weighting
   //    - This only seems to be used in
   //       - initial()
   //       - add_source(), but only for source_alpha()
   //       - output() to give a fiducial cell center
   //
   // ============================================= //

   double r3 = (rp*rp*rp + rp*rp*rm + rp*rm*rm + rm*rm*rm)/4.;
   double r2 = (rp*rp + rp*rm + rm*rm)/3.;
   double r = r3/r2;
   return( r );
}

double get_dV( double rp , double rm ){

   // ============================================= //
   //
   //  Calculates the volume of a cell
   //
   //  Inputs:
   //    - rp  - outer cell boundary
   //    - rm  - inner cell boundary
   //
   //  Returns:
   //    - dV  - cell volume
   //
   //  Side effects:
   //    None
   //
   //  Notes:
   //
   // ============================================= //

   // // this does the same thing
   // //   it's more consistent with the notation above
   // //   but less clear about what's happening physically
   // double dr  = rp-rm;
   // double r2    = (rp*rp+rm*rm+rp*rm)/3.;
   // return( 4.*M_PI*r2*dr );

   return( (4.*M_PI/3.)*(pow(rp,3.) - pow(rm,3.)) );

}
