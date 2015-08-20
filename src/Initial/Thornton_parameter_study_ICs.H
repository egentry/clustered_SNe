
#ifndef _THORNTON_PARAMETER_STUDY_ICS_H_
#define _THORNTON_PARAMETER_STUDY_ICS_H_

#include "initial_conditions.H"
#include "../structure.H"

class Thornton_Parameter_Study_ICs : public Initial_Conditions
{


public:
    Thornton_Parameter_Study_ICs(){E_blast=1e51;}

    int setICparams( struct domain * theDomain );
    void initial( double * prim , double r );
    int parse_command_line_args(    struct domain * theDomain , 
                                    int argc , 
                                    char * argv [] );

private:

    int setup_parameter_study( struct domain * theDomain );


    double mu;     // mean molecular weight -- this is the only time it's used
    double Gamma;  // adiabatic index

    double E_blast;        // [erg]

    int    completed_runs; // for starting a parameter study mid-way

    double background_density;      // [g cm^-3]
    double background_temperature;  // [K]
    double metallicity;

};

#endif