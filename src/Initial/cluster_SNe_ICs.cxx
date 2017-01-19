
#include <stdio.h>
#include <cmath>
#include <string>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <random>
#include <stdexcept>

#include "../misc.H" // get_mean_molecular_weight
#include "../structure.H"
#include "../constants.H"
#include "../blast.H"
#include "../mass_loss.H"
#include "initial_conditions.H"
#include "cluster_SNe_ICs.H"

const std::string Cluster_SNe_ICs::class_name = "cluster_SNe";

Cluster_SNe_ICs::Cluster_SNe_ICs(const double E_blast) 
    :   Initial_Conditions( class_name ),
        E_blast(E_blast)
{}

Cluster_SNe_ICs::Cluster_SNe_ICs(std::string overwrite_name) 
    :   Initial_Conditions( overwrite_name ),
        E_blast(1e51)
{}

int Cluster_SNe_ICs::setICparams( struct domain * theDomain ,
                                  const Mass_Loss * mass_loss)
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

    metallicity            = theDomain->metallicity;
    background_density     = theDomain->background_density;
    background_temperature = theDomain->background_temperature;

    mu = get_mean_molecular_weight( theDomain->metallicity );

    return 0;
}


void Cluster_SNe_ICs::initial( double * prim , double r )
{
    // initial background conditions
    const double rho = background_density;
    const double T_0 = background_temperature;  // [K]
    const double P   = (rho / (mu * m_proton)) * k_boltzmann * T_0;

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[VRR] = 0.0;
    prim[ZZZ] = metallicity;

}


int Cluster_SNe_ICs::setup_parameter_study( struct domain * theDomain )
{
    // ============================================= //
    //
    //  Sets up the initial conditions of each processor 
    //     to roughly match the parameter space explored by 
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
                                       metallicity_solar * pow(10, -3.0) };

    const int    n_background_densities   = 6;
    const double   background_densities[] = { m_proton*1.33e0,
                                              m_proton*1.33e+1,
                                              m_proton*1.33e+2,
                                              m_proton*1.33e-1,
                                              m_proton*1.33e-2,
                                              m_proton*1.33e-3 };

    const int    n_cluster_masses   = 5;
    const double   cluster_masses[] = { pow(10, 2),
                                        pow(10, 2.5),
                                        pow(10, 3),
                                        pow(10, 4),
                                        pow(10, 5) };

    if ( completed_runs >= (n_metallicities * n_background_densities 
                                            * n_cluster_masses) )
    {
        throw std::out_of_range("`completed_runs` exceeds possible values ");
    }

    int run=0;
    for( int i=0 ; i<n_metallicities ; ++i )
    {
        for( int j=0 ; j<n_background_densities ; ++j )
        {
            for( int k=0 ; k<n_cluster_masses ; ++k )
            {
                if( completed_runs == run )
                {
                    theDomain->metallicity            = metallicities[i];
                    theDomain->background_density     = background_densities[j];
                    theDomain->background_temperature = 167.8950;
                    theDomain->cluster_mass           = cluster_masses[k];
                }
                ++run;
            }
        }
    }





    // Assumes E_blast = 1e51

    double R = 300 * pc 
         * pow(E_blast / 1e51, 2./7) 
         * pow(theDomain->background_density / (1.33*m_proton), -.33)
         * pow(theDomain->cluster_mass / pow(10,2), .33);

    theDomain->theParList.rmax = R;
    theDomain->theParList.rmin = theDomain->theParList.rmax / 1e4;

    return 0;

}

int Cluster_SNe_ICs::parse_command_line_args (  struct domain * theDomain , 
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

    std::random_device rd;
    theDomain->seed = rd();
    if ( argc > 3 )
    {
        char *buf;
        theDomain->seed = strtol( argv[3] , &buf, 10);
        printf("seed = %u \n", theDomain->seed);
    }

    int error = 0;
    error = this->setup_parameter_study( theDomain );
    if ( error>0 ) return(error);

    return 0;
}


void Cluster_SNe_ICs::add_SNe( struct domain * theDomain,
                               const Mass_Loss * mass_loss )
{
    // maybe there's a better spot for this?
    theDomain->SNe = get_SNe(theDomain->cluster_mass,
                             theDomain->metallicity,
                             theDomain->seed,
                             mass_loss);
}


void Cluster_SNe_ICs::set_times( struct domain * theDomain )
{

    assert( std::is_sorted( theDomain->SNe.rbegin(),
                            theDomain->SNe.rend(),
                            sort_by_lifetime) );

    if ( theDomain->SNe.size() > 0 )
    {

        double t_first_SN = theDomain->SNe.back().lifetime;
        double t_last_SN  = theDomain->SNe.front().lifetime;

        theDomain->t      += t_first_SN;
        theDomain->t_init += t_first_SN;
        theDomain->t_fin   = 3 * t_last_SN;

    }
    else
    {
        // No SNe, so trust the input file for start and end times
    }

}

