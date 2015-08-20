
#ifndef _RESTART_ICS_H_
#define _RESTART_ICS_H_

#include <string>
#include "initial_conditions.H"
#include "../structure.H"

class Restart_ICs : public Initial_Conditions
{

public:
    Restart_ICs();
    ~Restart_ICs();

    int setICparams( struct domain * theDomain );
    void initial( double * prim , double r );
    int parse_command_line_args ( struct domain * theDomain , int argc , 
                                                              char * argv [] );
    void setup_grid( struct domain * theDomain );

private:

    int NL; // number of lines to be read in
    double * rr ;
    double * rho;
    double * Pp ;
    double * vr ;
    double * Z  ;
    double * A  ;

    std::string restart_filename;
    std::string restart_id;

    double time_restart;
    unsigned int checkpoints_finished;

    int countlines( std::string filename );

    int get_table( std::string filename );

    std::string get_restart_filename( const std::string partial_ID );

    void free_table();
};

#endif