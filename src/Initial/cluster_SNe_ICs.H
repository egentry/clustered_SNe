
#ifndef _CLUSTER_SNE_ICS_H_
#define _CLUSTER_SNE_ICS_H_

#include <string>

#include "initial_conditions.H"
#include "../structure.H"
#include "../mass_loss.H"

class Cluster_SNe_ICs : public Initial_Conditions
{

public:
    Cluster_SNe_ICs(const double E_blast = 1e51);
    Cluster_SNe_ICs(std::string overwrite_name);

    int virtual setICparams( struct domain * theDomain ,
                     const Mass_Loss * mass_loss );
    void virtual initial( double * prim , double r );
    int virtual parse_command_line_args (   struct domain * theDomain , 
                                    int argc , 
                                    char * argv [] );

    static const std::string class_name;
    
protected:
    void virtual add_SNe( struct domain * theDomain ,
                  const Mass_Loss * mass_loss );

    int virtual setup_parameter_study( struct domain * theDomain );

    double mu;
    double Gamma;  // adiabatic index

    double E_blast;        // [erg]

    int    completed_runs; // for starting a parameter study mid-way

    double background_density;      // [g cm^-3]
    double background_temperature;  // [K]
    double metallicity;             // mass fraction

private:


};

#endif