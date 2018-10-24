
#include <string>
#include <iostream>
#include <stdexcept>


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

extern "C" {
#include <grackle.h>
}

#include "constants.H" // defines physical constants
#include "cooling.H"
#include "structure.H" 

#include <boost/algorithm/string.hpp> 
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;   



Cooling * select_cooling( std::string cooling_type , const bool with_cooling )
{

    std::cout << "Using cooling type: " << cooling_type << std::endl;

    boost::algorithm::to_lower(cooling_type);

    std::string cooling_type_tmp;

    cooling_type_tmp = "equilibrium";
    boost::algorithm::to_lower(cooling_type_tmp);
    if ( cooling_type.compare(cooling_type_tmp) == 0 )
    {
        return new Equilibrium_Cooling (with_cooling);
    }

    throw std::invalid_argument("Cooling type didn't match known options");

}

//********************** Cooling (abstract base class) *********************//
Cooling::Cooling( const bool with_cooling , const std::string name)
    :   with_cooling(with_cooling),
        name(name)
{}


//********************** Equilibrium_Cooling *********************//
const std::string Equilibrium_Cooling::class_name = "equilibrium";


Equilibrium_Cooling::Equilibrium_Cooling( bool with_cooling ) 
    :   Cooling(with_cooling,
        class_name),
        grid_rank(1),
        field_size(1)
{

    // for ( int i=0 ; i<3 ; ++i )
    // {
    //     grid_start[i] = 0;
    //     grid_end[i] = field_size - 1;
    //     grid_dimension[i] = 1;
    // }

}

double Equilibrium_Cooling::calc_cooling( const double * prim , const double * cons , 
                                          const double dt, const bool cached )
{

    if( cached )
    {
        return dE_saved;
    }

    // declare fluid variable arrays (will need more for other chemistries)
    // gr_float *density, *energy, *x_velocity, *y_velocity, *z_velocity;
    // gr_float *metal_density;

    grackle_field_data my_fields;

    my_fields.grid_dimension = new int[3];
    my_fields.grid_start = new int[3];
    my_fields.grid_end = new int[3];

    my_fields.density         = new gr_float[field_size];
    my_fields.internal_energy = new gr_float[field_size];
    my_fields.x_velocity      = new gr_float[field_size];
    my_fields.y_velocity      = new gr_float[field_size];
    my_fields.z_velocity      = new gr_float[field_size];
    my_fields.metal_density   = new gr_float[field_size];

    const double density_initial = prim[RHO];
    const double energy_initial  = (1. / (grackle_data->Gamma - 1.)) * prim[PPP] / prim[RHO]
        / (cooling_units.velocity_units * cooling_units.velocity_units);     // internal energy PER UNIT MASS

    //copy old information into gr_float arrays
    for( int i=0 ; i<field_size ; ++i )
    {
        my_fields.density[i]          = density_initial / cooling_units.density_units;
        my_fields.internal_energy[i]  = energy_initial ;     // units above
        my_fields.x_velocity[i]       = prim[VRR] / cooling_units.velocity_units; // radial velocity
        my_fields.y_velocity[i]       = 0;
        my_fields.z_velocity[i]       = 0;
        my_fields.metal_density[i]    = prim[ZZZ] * my_fields.density[i];
    }
    for (int i = 0;i < 3;i++) {
        my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
        my_fields.grid_start[i] = 0;
        my_fields.grid_end[i] = field_size-1;
    }

    if (solve_chemistry(&cooling_units, &my_fields,
                                dt / cooling_units.time_units
                                ) == 0)
    {
        std::runtime_error("Error in solve_chemistry.");
    }

    // energy per unit volume
    const double dE = (my_fields.internal_energy[0] - energy_initial) * density_initial 
        * cooling_units.velocity_units * cooling_units.velocity_units;
    dE_saved = dE;

    delete my_fields.grid_dimension;
    delete my_fields.grid_start;
    delete my_fields.grid_end;
    delete my_fields.density;
    delete my_fields.internal_energy;
    delete my_fields.x_velocity;
    delete my_fields.y_velocity;
    delete my_fields.z_velocity;
    delete my_fields.metal_density;

    return dE;

};    

void Equilibrium_Cooling::setup_cooling( const struct domain * theDomain )
{
    if ( with_cooling == false) return;

    printf("setting up cooling \n");

    my_grackle_data = new chemistry_data;
    if (set_default_chemistry_parameters(my_grackle_data) == 0)
    {
        std::runtime_error("Error in set_default_chemistry_parameters");
    }



    char grackle_data_file[80] = "";
    strcat(grackle_data_file, GRACKLE_DIR);  // "GRACKLE_DIR" set by preprocessor
    strcat(grackle_data_file, "/input/CloudyData_UVB=HM2012.h5");
    // strcat(grackle_data_file, "/input/CloudyData_noUVB.h5");

    // check that file exists
    fs::path grackle_data_path(grackle_data_file);
    assert(fs::exists(grackle_data_path));
    printf("grackle data file: %s \n", grackle_data_file);

    // 0 = False,  1 = True
    // Many of these are the defaults. 
    // See: https://grackle.readthedocs.org/en/latest/Parameters.html#parameters

    grackle_data->use_grackle            = 1;  // turn on cooling
    grackle_data->with_radiative_cooling = 1;  // use cooling when updating chemistry
    grackle_data->primordial_chemistry   = 0;  // no chemistry network
    grackle_data->h2_on_dust             = 0;  // Don't use Omukai (2000)
    grackle_data->metal_cooling          = 1;  // allow metal cooling
    // grackle_data->cmb_temperature_floor  = 1;  // don't allow cooling below CMB temp.
    grackle_data->UVbackground           = 1;  // include UV background from grackle_data_file
    grackle_data->grackle_data_file      = grackle_data_file; // Haardt + Madau (2012)
    grackle_data->Gamma                  = theDomain->theParList.Adiabatic_Index;
    // All others are defaults...


    cooling_units.comoving_coordinates = 0;
    cooling_units.density_units        = m_proton; // need to ensure density field ~ 1
    cooling_units.time_units           = 1e8; // s 
    cooling_units.length_units         = 1e6 * cooling_units.time_units; // cm -- to get vel units of 10km/s
    cooling_units.velocity_units       = cooling_units.length_units / cooling_units.time_units;
    cooling_units.a_units              = 1.0;

    cooling_units.a_value = 1. / (1. + theDomain->theParList.Cooling_Redshift);

    if (initialize_chemistry_data(&cooling_units) == 0)
    {
        std::runtime_error("Error in initialize_chemistry_data");
    }


    printf("use_grackle: %d \n", grackle_data->use_grackle);
    grackle_verbose=1;


};
