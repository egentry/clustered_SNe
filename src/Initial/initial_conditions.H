
#ifndef _INITIAL_CONDITIONS_H_
#define _INITIAL_CONDITIONS_H_

#include <vector>
#include <string>

#include "../structure.H"

#include "initial_conditions.H"


class Initial_Conditions 
{

public:
    virtual int  setICparams( struct domain * theDomain ){return 0;}
    virtual int  parse_command_line_args(   struct domain * theDomain , 
                                            int argc , 
                                            char * argv [] ){return 0;}

    virtual void setup_grid( struct domain * theDomain );
    virtual void possibly_extend_grid( struct domain * theDomain);
    virtual void extend_grid( struct domain * theDomain ,
                              const double extend_fraction=0.25);

    virtual bool trust_LogZoning_flag(){return true;}

    virtual void set_times( struct domain * theDomain );

private:
    virtual void initial( double * prim , double r ) = 0;
    virtual void setup_cells( struct domain * theDomain );

};

Initial_Conditions * select_initial_conditions( std::string IC_name );

#endif

