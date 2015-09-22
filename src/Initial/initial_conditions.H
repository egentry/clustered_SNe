
#ifndef _INITIAL_CONDITIONS_H_
#define _INITIAL_CONDITIONS_H_

#include <vector>
#include <string>

#include "../structure.H"

#include "../mass_loss.H" // Mass_Loss


class Initial_Conditions 
{

public:
    Initial_Conditions( const std::string name );

    virtual int  setICparams( struct domain * theDomain ,
                              const Mass_Loss * mass_loss );
    virtual int  parse_command_line_args( struct domain * theDomain , 
                                          int argc , 
                                          char * argv [] );

    virtual void setup_grid( struct domain * theDomain );
    virtual void possibly_extend_grid( struct domain * theDomain);
    virtual void extend_grid( struct domain * theDomain ,
                              const double extend_fraction=0.25);

    virtual bool trust_LogZoning_flag() const ;

    const std::string name;

protected:
    virtual void add_SNe( struct domain * theDomain ,
                          const Mass_Loss * mass_loss );

    virtual void set_times( struct domain * theDomain );
    virtual void set_output_prefix( struct domain * theDomain );

private:
    virtual void initial( double * prim , double r ) = 0;
    virtual void setup_cells( struct domain * theDomain );

};

Initial_Conditions * select_initial_conditions( std::string IC_name );

#endif

