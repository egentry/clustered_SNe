
#ifndef _COOLING_H_
#define _COOLING_H_

#include <string>

extern "C" {
#include <grackle.h>
}

class Cooling
{


public:
    Cooling( const bool with_cooling , const std::string name );
    virtual double calc_cooling( const double * prim , const double * cons , 
                                 const double dt, const bool cached ) = 0;
    virtual void setup_cooling(  const struct domain * theDomain ) = 0;  

    const bool with_cooling;
    const std::string name;

protected:

private:


};

class Equilibrium_Cooling : public Cooling
{
public:
    
    Equilibrium_Cooling( bool with_cooling );

    double calc_cooling( const double * prim , const double * cons , 
                         const double dt, const bool cached );

    void setup_cooling( const struct domain * theDomain );

    static const std::string class_name;

protected:

private:
    chemistry_data *my_grackle_data;
    code_units cooling_units;

    const int grid_rank;
    const int field_size;

    // double a_value;

    // these can't be constantified due to grackle
    // int grid_start[3];
    // int grid_end[3];
    // int grid_dimension[3];

    double dE_saved;

};

Cooling * select_cooling( std::string cooling_type ,
                          const bool with_cooling );
          

#endif 
