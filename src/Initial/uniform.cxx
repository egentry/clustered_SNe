
#include "../structure.H"

int setICparams( struct domain * theDomain ){
	return(0);
}

void initial( double * prim , double r ){
   prim[RHO] = 1.0;
   prim[PPP] = 1.0;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;
}

int parse_command_line_args ( struct domain * theDomain , int argc , char * argv [] )
{
      return(0);
}