
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <assert.h>

#include <algorithm>
#include <vector>

#include <random>

#include "../structure.H"
#include "../constants.H"
#include "../blast.H"
#include "../mass_loss.H"
#include "initial_conditions.H"
#include "Thornton_parameter_study_ICs.H"


const std::string Thornton_Parameter_Study_ICs::class_name = "disappear";

Thornton_Parameter_Study_ICs::Thornton_Parameter_Study_ICs(const double E_blast) 
    :   Initial_Conditions( class_name ),
        E_blast(E_blast)
{}

int Thornton_Parameter_Study_ICs::setICparams( struct domain * theDomain ,
                                               const Mass_Loss * mass_loss )
{

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

    this->set_output_prefix( theDomain );
    this->add_SNe( theDomain , mass_loss );
    this->set_times( theDomain );

    Gamma = theDomain->theParList.Adiabatic_Index;

    background_density     = theDomain->background_density;
    background_temperature = theDomain->background_temperature;

    const double Y = .23; // helium fraction
    const double Z = theDomain->metallicity; // metals fraction
    const double X = 1 - Y - Z; // hydrogen mass fraction

    mu = 1. / (2*X + .75*Y + .5*Z); // mean molecular weight

    return 0;
}


void Thornton_Parameter_Study_ICs::initial( double * prim , double r )
{
    double rho,Pp;

    // initial background conditions
    rho = background_density;
    const double T_0 = background_temperature;  // [K]
    Pp  = (rho / (mu * m_proton)) * k_boltzmann * T_0;

    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[VRR] = 0.0;
    prim[ZZZ] = metallicity;

}

int Thornton_Parameter_Study_ICs::setup_parameter_study( struct domain * theDomain )
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
        std::cerr << "ERROR: enter run does not correspond to a Thornton IC" 
                  << std::endl;
        return 1;
    }

    int run=0;
    for( int i=0 ; i<n_metallicities ; ++i )
    {
        for( int j=0 ; j<n_background_densities ; ++j )
        {
            if( completed_runs == run )
            {
                theDomain->metallicity            = metallicities[i];
                metallicity                       = metallicities[i];
                theDomain->background_density     = background_densities[j];
                theDomain->background_temperature = 1e4;
            }
            ++run;
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

    return 0;

}

int Thornton_Parameter_Study_ICs::parse_command_line_args( struct domain * theDomain , 
                                                           int argc , 
                                                           char * argv [] )
{
    // ============================================= //
    //
    //  Parses the command line args, as necessary into initial condition info
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
    if ( argc > 2 )
    {
        char *buf;
        completed_runs = strtol( argv[2] , &buf, 10);
        printf("completed_runs = %d \n", completed_runs);
    }

    int error = 0;
    error = this->setup_parameter_study( theDomain );
    if ( error>0 ) return(error);

    return 0;
}

void Thornton_Parameter_Study_ICs::add_SNe( struct domain * theDomain,
                                            const Mass_Loss * mass_loss )
{
    // set only one supernova
    supernova tmp;
    tmp.lifetime = 0.0;
    tmp.mass = 0;
    tmp.mass_ejecta   = 3 * M_sun;
    tmp.mass_ejecta_Z = tmp.mass_ejecta * theDomain->metallicity;
    tmp.mass_winds = 0;
    theDomain->SNe.push_back(tmp);
}