
#ifndef _SEDOV_ICS_H_
#define _SEDOV_ICS_H_

#include <string>

#include "initial_conditions.H"
#include "../structure.H"
#include "../mass_loss.H"


class Sedov_ICs : public Initial_Conditions
{


public:

    Sedov_ICs();

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

private:

    double mu;
    double Gamma;  // adiabatic index

    double E_blast;        // [erg]

    double background_density;      // [g cm^-3]
    double background_temperature;  // [K]
    double metallicity;             // mass fraction

};

#endif