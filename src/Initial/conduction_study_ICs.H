
#ifndef _CONDUCTION_STUDY_ICS_H_
#define _CONDUCTION_STUDY_ICS_H_

#include <string>

#include "initial_conditions.H"
#include "cluster_SNe_ICs.H"

#include "../structure.H"
#include "../mass_loss.H"

class Conduction_Study_ICs : public Cluster_SNe_ICs
{

public:
    Conduction_Study_ICs();

    int virtual parse_command_line_args (   struct domain * theDomain , 
                                    int argc , 
                                    char * argv [] );

    static const std::string class_name;
    
protected:


private:


};

#endif