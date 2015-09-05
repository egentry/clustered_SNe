
#include <iostream>

#include "structure.H"

#include "blast.H"
#include "readpar.H" // read_par_file()
#include "domain.H" // check_dt, possiblyOutput, setupDomain, freeDomain
#include "misc.H" // getmindt, set_wcell
#include "mass_loss.H" // Mass_Loss, add_mass_loss
#include "profiler.H" // start_clock, generate_log
#include "timestep.H"
#include "Output/ascii.H" // overview
#include "Initial/initial_conditions.H" // select_initial_conditions

int main( int argc , char * argv[] )
{

    int error;

    struct domain theDomain = {0};

    error = read_par_file( &theDomain , argc , argv );
    if( error==1 ) 
    {
        std::cerr << "Error in read_par_file" << std::endl;
        return 0;
    }

    Initial_Conditions * ICs;
    ICs = select_initial_conditions(theDomain.theParList.ICs);

    Mass_Loss * mass_loss;
    mass_loss = select_mass_loss(theDomain.theParList.mass_loss);

    error = ICs->parse_command_line_args( &theDomain , argc , argv );
    if( error==1 ) 
    {
        std::cerr << "Error in parse_command_line_args" << std::endl;
        return 0;
    }

    ICs->add_SNe( &theDomain , mass_loss );

    error = setupDomain( &theDomain , ICs , mass_loss );
    if( error==1 ) 
    {
        std::cerr << "Error in setupDomain" << std::endl;
        return 0;
    }

    ICs->setup_grid( &theDomain );   

    error = overview( &theDomain );
    if( error==1)
    {
        std::cerr << "Error in overview" << std::endl;
        return 0;
    }

    start_clock( &theDomain ); 

    while( !(theDomain.final_step) )
    {
        add_blasts( &theDomain );
        set_wcell( &theDomain );
        double dt = getmindt( &theDomain );
        check_dt( &theDomain , &dt );
        possiblyOutput( &theDomain , 0 );
        timestep( &theDomain , dt , ICs );
        mass_loss->add_mass_loss( &theDomain , dt );

    }

    possiblyOutput( &theDomain , 1 );

    generate_log( &theDomain );
    freeDomain( &theDomain );

    return 0;

}