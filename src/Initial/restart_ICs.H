
#ifndef _RESTART_ICS_H_
#define _RESTART_ICS_H_

#include <string>
#include "initial_conditions.H"
#include "../structure.H"
#include "../mass_loss.H"

class Restart_ICs : public Initial_Conditions
{

public:
    Restart_ICs();
    ~Restart_ICs();

    int setICparams( struct domain * theDomain , 
                     const Mass_Loss * mass_loss);
    void initial( double * prim , double r );
    int parse_command_line_args ( struct domain * theDomain , int argc , 
                                                              char * argv [] );
    void setup_grid( struct domain * theDomain );
    bool trust_LogZoning_flag() const ;

    static const std::string class_name;

protected:
    void add_SNe( struct domain * theDomain ,
                  const Mass_Loss * mass_loss );

    void set_times( struct domain * theDomain );
    void set_output_prefix( struct domain * theDomain );

private:

    int NL; // number of lines to be read in
    double * rr ;
    double * rho;
    double * Pp ;
    double * vr ;
    double * Z  ;

    // keller-specific variables
    int * multiphase ;
    double * x_cold ;
    double * y_cold ;
    double * z_cold ;
    double * x_hot ;
    double * y_hot ;
    double * z_hot ;


    std::string restart_filename;
    std::string restart_id;

    int N_chk;
    double delta_time;
    double CFL;
    double Cooling_Redshift;


    double time_restart;
    unsigned int last_finished_checkpoint;

    int get_table( const std::string filename );
    std::string get_restart_filename( const std::string partial_ID );
    std::string filename_to_prefix( const std::string filename ) const;

    void free_table();
};

#endif