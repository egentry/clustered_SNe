
#include <iostream>
#include <string>
#include <random>


#include "initial_conditions.H"
#include "cluster_SNe_ICs.H"
#include "conduction_study_ICs.H"


const std::string Conduction_Study_ICs::class_name = "conduction";

Conduction_Study_ICs::Conduction_Study_ICs() 
    :   Cluster_SNe_ICs( class_name )
{}

int Conduction_Study_ICs::parse_command_line_args (  struct domain * theDomain , 
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

    if ( argc > 4 )
    {
        theDomain->theParList.H_0 = stod(std::string(argv[4]));
        printf("H_0 = %lf \n", theDomain->theParList.H_0);
    }
    else
    {
        std::cerr << "Missing argument " << 3 << std::endl;
        return 1;
    }

    if ( argc > 5 )
    {
        theDomain->theParList.H_1 = stod(std::string(argv[5]));
        printf("H_1 = %lf \n", theDomain->theParList.H_1);
    }
    else
    {
        std::cerr << "Missing argument " << 4 << std::endl;
        return 1;
    }

    int error = 0;
    error = this->setup_parameter_study( theDomain );
    if ( error>0 ) return(error);

    return 0;
}