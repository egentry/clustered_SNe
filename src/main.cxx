
#include <iostream>

#include "structure.H"

#include "blast.H"
#include "cooling.H" // Cooling
#include "readpar.H" // read_par_file()
#include "domain.H" // check_dt, possiblyOutput, setupDomain, freeDomain
#include "misc.H" // getmindt, set_wcell
#include "mass_loss.H" // Mass_Loss, add_mass_loss
#include "profiler.H" // start_clock, generate_log
#include "timestep.H"
#include "Output/ascii.H" // overview
#include "Initial/initial_conditions.H" // select_initial_conditions
#include <cblas.h>

int main( int argc , char * argv[] )
{

    #ifdef _OPENMP
        printf("Using OpenMP\n");
    #endif

    int error;

    goto_set_num_threads( NThread );
    openblas_set_num_threads( NThread );

    struct domain theDomain = {0};

    error = read_par_file( &theDomain , argc , argv );
    if( error==1 ) 
    {
        std::cerr << "Error in read_par_file" << std::endl;
        return 0;
    }

    Initial_Conditions * ICs;
    ICs = select_initial_conditions(theDomain.theParList.ICs);


    error = ICs->parse_command_line_args( &theDomain , argc , argv );
    if( error==1 ) 
    {
        std::cerr << "Error in parse_command_line_args" << std::endl;
        return 0;
    }

    // initial conditions might change the cooling and mass_loss types


    Mass_Loss * mass_loss;
    mass_loss = select_mass_loss(theDomain.theParList.mass_loss);

    error = setupDomain( &theDomain , ICs , mass_loss );
    if( error==1 ) 
    {
        std::cerr << "Error in setupDomain" << std::endl;
        return 0;
    }

    ICs->setup_grid( &theDomain );   


    // cooling might depend on the actual grid
    Cooling * cooling;
    cooling = select_cooling(theDomain.theParList.cooling_type , 
                             theDomain.theParList.With_Cooling );
    cooling->setup_cooling( &theDomain ); // put in constructor? 

    overviews( &theDomain , mass_loss , cooling );

    start_clock( &theDomain ); 

    while( !(theDomain.final_step) )
    {
        add_blasts( &theDomain );
        set_wcell( &theDomain );
        double dt = getmindt( &theDomain );
        check_dt( &theDomain , &dt );
        possiblyOutput( &theDomain , 0 );
        timestep( &theDomain , dt , ICs , cooling );
        mass_loss->add_mass_loss( &theDomain , dt );

    }

    possiblyOutput( &theDomain , 1 );

    generate_log( &theDomain );
    freeDomain( &theDomain );

    return 0;

}