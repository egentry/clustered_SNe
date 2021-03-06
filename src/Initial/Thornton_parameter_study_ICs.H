
#ifndef _THORNTON_PARAMETER_STUDY_ICS_H_
#define _THORNTON_PARAMETER_STUDY_ICS_H_

#include <string>

#include "initial_conditions.H"
#include "../structure.H"
#include "../mass_loss.H"


class Thornton_Parameter_Study_ICs : public Initial_Conditions
{


public:

    Thornton_Parameter_Study_ICs(const double E_blast = 1e51);

    int setICparams( struct domain * theDomain ,
                     const Mass_Loss * mass_loss );
    void initial( double * prim , double r );
    int parse_command_line_args(    struct domain * theDomain , 
                                    int argc , 
                                    char * argv [] );

    static const std::string class_name;

protected:
    void add_SNe( struct domain * theDomain ,
                  const Mass_Loss * mass_loss );

    int setup_parameter_study( struct domain * theDomain );

private:

    double mu;
    double Gamma;  // adiabatic index

    double E_blast;        // [erg]

    int    completed_runs; // for starting a parameter study mid-way

    double background_density;      // [g cm^-3]
    double background_temperature;  // [K]
    double metallicity;             // mass fraction

};

#endif