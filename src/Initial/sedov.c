
#include <stdio.h>
#include <math.h>

#include "../structure.h"
#include "../constants.h"


// ======== STATIC VARIABLES SPECIFIC TO THESE INITIAL CONDITIONS ========= //

 
static double mu;     // mean molecular weight -- this is the only time it's used
static double Gamma;  // adiabatic index

static double E_blast = 1e51;        // [erg]
static double M_blast = 3 * M_sun;   // [g]
static double V_blast;  // volume of the initial blast (1 cell)
static double R_blast;  // outermost boundary of initial blast (1 cell)

static int    completed_runs; // for starting a parameter study mid-way

static double background_density;      // [g cm^-3]
static double background_temperature;  // [K]

// ======== Setup functions ========= //


double get_dV( double , double );
int setup_parameter_study( struct domain * theDomain );
int setICparams( struct domain * theDomain ){

   // ============================================= //
   //
   //  Doesn't set the actual initial conditions;
   //    simply sets the parameters used for constructing the initial conditions
   //    (which will be stored as static variables within the scope of this file)
   //    since we won't have access to the entire domain
   //    when we are later calling the initial() function
   //
   //  Inputs:
   //     - theDomain    - the standard domain struct used throughout
   //
   //  Outputs:
   //     - error        - 0 for successful execution
   //                    - 1 for failure (should cause this run to quietly stop)
   //
   //  Side effects:
   //     - overwrites the static variables for these initial conditions
   //
   //  Notes:
   //     - By this point the grid spacing has already been set,
   //       but not the fluid variables at those grid points
   //
   // ============================================= //



   Gamma = theDomain->theParList.Adiabatic_Index;

   int i;
   const int n_blast = 2;
   V_blast = 0.;
   R_blast = 0.;
   for( i=1 ; i<n_blast ; ++i )
   {
      struct cell * c = &(theDomain->theCells[i]);
      const double rp = c->riph;
      const double rm = rp - c->dr; 
      V_blast += get_dV( rp , rm );
      R_blast = rp;
   }

   background_density     = theDomain->background_density;
   background_temperature = theDomain->background_temperature;

   const double Y = .23; // helium fraction
   const double Z = theDomain->metallicity; // metals fraction
   const double X = 1 - Y - Z; // hydrogen mass fraction

   mu = 1. / (2*X + .75*Y + .5*Z); // mean molecular weight

   return(0);
}

void initial( double * prim , double r ){
   double rho,Pp;
   double X;

   if( r <= R_blast ){


      // rho = M_blast / V_blast;
      rho = background_density;

      Pp  = (Gamma-1) * E_blast / V_blast;    
      X = 1.0;
   }else{
      // initial background conditions
      rho = background_density;
      const double T_0 = background_temperature;  // [K]
      Pp  = (rho / (mu * m_proton)) * k_boltzmann * T_0;

      X = 0.0;
   }
   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = 0.0;
   prim[XXX] = X;
   prim[AAA] = 0.0;
}



int setup_parameter_study( struct domain * theDomain )
{
   // ============================================= //
   //
   //  Sets up the initial conditions of each processor 
   //     to match the parameter space explored by 
   //     Thornton et al. (1998), ApJ, 500, 95
   //
   //  Inputs:
   //     - theDomain    - the standard domain struct used throughout
   //
   //  Outputs:
   //     - error        - 0 for successful execution
   //                    - 1 for failure (should cause this run to quietly stop)
   //
   //  Side effects:
   //     - overwrites:
   //        - theDomain->metallicity
   //        - theDomain->background_density
   //        - theDomain->background_temperature
   //
   //  Notes:
   //     - Doesn't actually scale t_end, r_max appropriately
   //        - t_end could be changed here, but r_max is already set
   //     - If you need to one/some runs over,
   //       you'll need to set completed_runs appropriately
   //        - later could make this a command line argument
   //
   //
   // ============================================= //


   const double metallicity_solar = .02;
   const int    n_metallicities   = 7;
   const double   metallicities[] = { metallicity_solar * pow(10,  0.0),
                                      metallicity_solar * pow(10,  0.5),
                                      metallicity_solar * pow(10, -0.5),
                                      metallicity_solar * pow(10, -1.0),
                                      metallicity_solar * pow(10, -1.5),
                                      metallicity_solar * pow(10, -2.0),
                                      metallicity_solar * pow(10, -3.0)};
   const int    n_background_densities   = 7;
   const double   background_densities[] = { m_proton*1.33e0,
                                             m_proton*1.33e+1,
                                             m_proton*1.33e+2,
                                             m_proton*1.33e+3,
                                             m_proton*1.33e-1,
                                             m_proton*1.33e-2,
                                             m_proton*1.33e-3};


   if( completed_runs >= (n_metallicities * n_background_densities) )
   {
      return(1);
   }

   int i,j,k;
   k=0;
   for( i=0 ; i<n_metallicities ; ++i )
   {
      for( j=0 ; j<n_background_densities ; ++j)
      {
         if( completed_runs == k )
         {
            theDomain->metallicity            = metallicities[i];
            theDomain->background_density     = background_densities[j];
            theDomain->background_temperature = 1e4;
         }
         ++k;
      }
   }


   // Set the scale of the simulation -- see Thornton et al. (1998) Eqs 14-35
   // Assumes E_blast = 1e51
   double R_thornton;
   if( log10(theDomain->metallicity / metallicity_solar) > -2  )
   {
      // Equation 20
      R_thornton = 49.3 * pc 
         * pow(E_blast / 1e51, 2./7) 
         * pow(theDomain->background_density / m_proton, -.42)
         * pow(theDomain->metallicity / metallicity_solar, -.1);
   }
   else
   {
      // Equation 31
      R_thornton = 78.1 * pc 
         * pow(E_blast / 1e51, 2./7) 
         * pow(theDomain->background_density / m_proton, -.42);
   }

   theDomain->theParList.rmax = 4 * R_thornton;
   theDomain->theParList.rmin = theDomain->theParList.rmax / 1e4;

   printf("R_min = %le \n", theDomain->theParList.rmin);
   printf("R_max = %le \n", theDomain->theParList.rmax);

   // Sets the end time appropriately,
   // having done a 2d power law fit to t_f(n_0, Z)
   // using the results of Thornton (Table 3)
   double t_f = 5.52e5 * yr 
         * pow(theDomain->background_density / m_proton,   -.53)
         * pow(theDomain->metallicity / metallicity_solar, -.16);

   theDomain->theParList.t_max = 5 * t_f;
   printf("t_max = %le \n", t_f);

   return(0);

}

int parse_command_line_args ( struct domain * theDomain , int argc , char * argv [] )
{
   // ============================================= //
   //
   //  Parses the command line args, as necessary into initial condition information
   //
   //  Inputs:
   //     - theDomain    - the standard domain struct used throughout
   //     - argc         - count of command line args
   //     - argv         - array of string pointers for each command line arg
   //
   //  Outputs:
   //     - error        - 0 for successful execution
   //                    - 1 for failure (should cause this run to quietly stop)
   //
   //  Side effects:
   //     - overwrites:
   //        - completed_runs (file-scope static variable)
   //        - 
   //
   //  Notes:
   //     - Could set r_max, t_end here if desired
   //       (would be useful for scaling domain to scale from inputs)
   //
   // ============================================= //

   completed_runs = 0;
   if ( argc > 1 )
   {
      char *buf;
      completed_runs = strtol( argv[1] , &buf, 10);
      printf("completed_runs = %d \n", completed_runs);
   }

   int error = 0;
   error = setup_parameter_study( theDomain );
   if ( error==1 ) return(error);

   return(0);
}