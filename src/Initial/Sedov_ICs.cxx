
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <assert.h>
#include <vector>
#include <string>

#include "../misc.H" // get_mean_molecular_weight
#include "../structure.H"
#include "../constants.H"
#include "../blast.H"
#include "../mass_loss.H"
#include "initial_conditions.H"
#include "Sedov_ICs.H"


const std::string Sedov_ICs::class_name = "Sedov";

Sedov_ICs::Sedov_ICs() 
    :   Initial_Conditions( class_name )
{}

int Sedov_ICs::setICparams( struct domain * theDomain ,
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

    mu = get_mean_molecular_weight( theDomain-> metallicity );

    return 0;
}


void Sedov_ICs::initial( double * prim , double r )
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


int Sedov_ICs::parse_command_line_args( struct domain * theDomain , 
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
    //     - sets:
    //        - background_density
    //        - metallicity
    //        - background temperature
    //        - E_blast
    //
    //  Notes:
    //     - argv[1] is reserved for the parameter file
    //
    // ============================================= //


    theDomain->background_density = m_proton;
    if ( argc > 2 )
    {
        theDomain->background_density = std::stod(std::string( argv[2] ));
        printf("background_density = %e \n", theDomain->background_density);
    }


    theDomain->metallicity = 0.02; // solar
    if ( argc > 3 )
    {
        theDomain->metallicity = std::stod(std::string( argv[3] ));
        printf("metallicity = %f \n", theDomain->metallicity);
    }


    theDomain->background_temperature = 1e4;
    if ( argc > 4 )
    {
        theDomain->background_temperature = std::stod(std::string( argv[4] ));
        printf("background_temperature = %e \n", theDomain->background_temperature);
    }


    E_blast = 1e51;
    if ( argc > 5 )
    {
        E_blast = std::stod(std::string( argv[5] ));
        printf("E_blast = %e \n", E_blast);
    }

    if (theDomain->theParList.With_Cooling == true)
    {
        std::cerr << "Warning: Overwriting With_Cooling from true to false" << std::endl;   
        // theDomain->theParList.With_Cooling = false;
    }

    return 0;
}

void Sedov_ICs::add_SNe( struct domain * theDomain,
                                            const Mass_Loss * mass_loss )
{
    // set only one supernova
    supernova tmp;
    tmp.lifetime = 0.0;
    tmp.mass = 0.0;
    tmp.mass_ejecta   = 0.0;
    tmp.mass_ejecta_Z = 0.0;
    tmp.mass_winds = 0.0;
    theDomain->SNe.push_back(tmp);
}