
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <grackle.h>

#include "constants.h" // defines physical constants
#include "structure.h" 

double calc_cooling( double * prim ,  double * cons , double metallicity , double dt , code_units my_units )
{

    int i;

    grackle_verbose=1;

    double z_redshift = 0.;
    double a_value = 1. / (1. + z_redshift);


    int field_size = 1; // we'll be removing guard cells + zone boundaries
    int grid_rank  = 1;
    int grid_dimension[3], grid_start[3], grid_end[3];
    for (i=0; i<3; ++i)
    {
        // set defaults -- for 1d only change element [0] later
        grid_dimension[i]   = 1;
        grid_start[i]       = 0;
        grid_end[i]         = 0;
    }
    // for the dimension we're using, set guard cell information + size
    grid_dimension[0]   = field_size;
    grid_start[0]       = 0;
    grid_end[0]         = (field_size - 1);


    // declare fluid variable arrays (will need more for other chemistries)
    gr_float *density, *energy, *x_velocity, *y_velocity, *z_velocity;
    gr_float *metal_density;

    density         = malloc(field_size * sizeof(gr_float));
    energy          = malloc(field_size * sizeof(gr_float));
    x_velocity      = malloc(field_size * sizeof(gr_float));
    y_velocity      = malloc(field_size * sizeof(gr_float));
    z_velocity      = malloc(field_size * sizeof(gr_float));
    metal_density   = malloc(field_size * sizeof(gr_float));

    double density_initial = prim[RHO];
    double energy_initial  = (1. / (grackle_data.Gamma - 1.)) * prim[PPP] / prim[RHO];

    //copy old information into gr_float arrays
    for(i=0; i<field_size; ++i)
    {
        density[i]          = density_initial / my_units.density_units;
        energy[i]           = energy_initial;     // internal energy PER UNIT MASS
        x_velocity[i]       = prim[VRR];          // radial velocity
        y_velocity[i]       = 0;
        z_velocity[i]       = 0;
        metal_density[i]    = metallicity * density[i];
    }


    // setup_cooling();

    // printf("done setting IC's \n");


    if (solve_chemistry_table(&my_units,
                                a_value, dt,
                                grid_rank, grid_dimension,
                                grid_start, grid_end,
                                density, energy,
                                x_velocity, y_velocity, z_velocity,
                                metal_density) == 0)
    {
        fprintf(stderr, "Error in solve_chemistry.\n");
    }

    // printf("done setting solving chemistry \n");

    double dE = (energy[0] - energy_initial) * density_initial; // energy per unit volume

    // printf("use_grackle: %d \n", grackle_data.use_grackle);
    free(density);
    free(energy);
    free(x_velocity);
    free(y_velocity);
    free(z_velocity);
    free(metal_density);
    return dE;

};    

code_units setup_cooling( struct domain * theDomain )
{

    printf("setting up cooling \n");


    if (set_default_chemistry_parameters() == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    }


    char grackle_data_file[80] = "";
    strcat(grackle_data_file, getenv("HOME"));
    strcat(grackle_data_file, "/local/grackle/input/CloudyData_UVB=HM2012.h5"); 

    printf("grackle data file: %s \n", grackle_data_file);

    // 0 = False,  1 = True
    // Many of these are the defaults. See: https://grackle.readthedocs.org/en/latest/Parameters.html#parameters

    grackle_data.use_grackle            = 1;  // turn on cooling
    grackle_data.with_radiative_cooling = 1;  // use cooling when updating chemistry
    grackle_data.primordial_chemistry   = 0;  // no chemistry network
    grackle_data.h2_on_dust             = 0;  // Don't use Omukai (2000)
    grackle_data.metal_cooling          = 1;  // allow metal cooling
    // grackle_data.cmb_temperature_floor  = 1;  // don't allow cooling below CMB temp.
    grackle_data.UVbackground           = 1;  // include UV background from grackle_data_file
    grackle_data.grackle_data_file      = grackle_data_file; // Haardt + Madau (2012)
    grackle_data.Gamma                  = theDomain->theParList.Adiabatic_Index;
    // All others are defaults...

    code_units my_units;
    my_units.comoving_coordinates = 0;
    my_units.density_units        = m_proton; // need to ensure density field ~ 1
    my_units.length_units         = 1.0;
    my_units.time_units           = 1.0;
    my_units.velocity_units       = my_units.length_units / my_units.time_units;
    my_units.a_units              = 1.0;
    double z_redshift = 0.;
    double a_value = 1. / (1. + z_redshift);

    if (initialize_chemistry_data(&my_units, a_value) == 0)
    {
        fprintf(stderr, "Error in initialize_chemistry_data.\n");
    }


    printf("use_grackle: %d \n", grackle_data.use_grackle);

    return my_units;
};