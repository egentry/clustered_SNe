#include <stdio.h>
#include <stdbool.h>
#include <math.h>

// #include <unistd.h>
#include <string.h>
#include <stdlib.h>

#include <grackle.h>

#include <constants.h> // defines physical constants
#include <fluid.adiabatic.h>
#include <cooling.h>
#include <ICs.h>
#include <grid.h>




// MOVE ALL THESE INTO AN INIT FILE


// total zone number variables must be ODD
int    zones       = 5*1001; // including guard cells (remember to add 4 guard cells)
int    zones_inner = 5*300; // not used for sedov ICs

char   IC[80]      = "sedov";

const bool   with_cooling = false;
const bool   with_gravity = false;


int main()
{

    int timesteps, n_print;

    timesteps   = 5*2000000; 
    n_print     = timesteps / 100; // print every n_steps
    // timesteps   = 10000; 
    // n_print     = 100; // print every n_steps

    const double metallicity = grackle_data.SolarMetalFractionByMass;

    const double a_visc           = 3;     // sets artificial viscosity strength
    const double eta              = 0.2;   // CFL number
          double eta_var          = 1.;    // Variable CFL knob
          bool   eta_var_continue = true;  // keep iterating to lower eta_var's?

    double       E_int, E_kin, E_grav, E_tot, Momentum; 

    double *U_old, *U_new;
    double *R_old, *R_new;
    double *V_old, *V_new;
    double *T_old, *T_new;
    double *E_old, *E_new, *dE;
    double *P_old, *P_new;
    double *Q_old, *Q_new;
    double *C_ad;
    double *Mass, *M_int;
    U_old = malloc(zones * sizeof(double));
    U_new = malloc(zones * sizeof(double));
    R_old = malloc(zones * sizeof(double));
    R_new = malloc(zones * sizeof(double));
    V_old = malloc(zones * sizeof(double));
    V_new = malloc(zones * sizeof(double));
    T_old = malloc(zones * sizeof(double));
    T_new = malloc(zones * sizeof(double));
    E_old = malloc(zones * sizeof(double));
    E_new = malloc(zones * sizeof(double));
    dE    = malloc(zones * sizeof(double));
    P_old = malloc(zones * sizeof(double));
    P_new = malloc(zones * sizeof(double));
    Q_old = malloc(zones * sizeof(double));
    Q_new = malloc(zones * sizeof(double));
    C_ad  = malloc(zones * sizeof(double));
    Mass  = malloc(zones * sizeof(double));
    M_int = malloc(zones * sizeof(double));

    int          i,j,k;  // loop variables
    int          k_start = 0; // start at timestep 0, unless continuing a run

    double       time_current = 0;

    // declare various temporary parameters
    double       time_delta_tmp, time_delta;
    double       accel_thermo, accel_grav;



    code_units my_units;
    if (with_cooling == true)
    {
        my_units = setup_cooling();
    }


    // ========== PREPARE OUTPUT FILES ============ //

    FILE *info_file_pointer, *velos_file_pointer; // outputs
    char   info_filename[80]   = "info";
    char  velos_filename[80]   = "velos";
    char    mod_filename[80]   = "modfile"; // used elsewhere


    if (with_cooling == true)
    {
        strcat(velos_filename, "_with_cooling");
        strcat( info_filename, "_with_cooling");
        strcat(  mod_filename, "_with_cooling");
    }

    printf("velos_filename: %s \n", velos_filename);


    // ========== MAKE INITIAL CONDITIONS ============ //

    if (strcmp(IC, "sedov")==0)
    {
        time_current = sedov(zones, U_old, R_old, V_old, T_old, 
            E_old, P_old, Q_old, C_ad); 
    }
    else if (strcmp(IC, "restart")==0)
    {
        time_current = restart(&zones, &k_start, U_old, R_old, V_old, T_old, 
            E_old, P_old, Q_old, C_ad, mod_filename); 
        printf("k_start = %d \n", k_start);
    }
    else
    {
        // INITIALIZE DATA FOR A NEW RUN
        new_blast(zones, zones_inner, U_old, R_old, V_old, T_old,
            E_old, P_old, Q_old, C_ad);
    }

    create_Mass_array(Mass,   R_old, V_old, zones);
    create_M_int_array(M_int, Mass, zones);


    // interpolate values as necessary
    for(i=1; i<zones-1; ++i)
    {
        if(i%2==0)
        {
            // ZONES
            R_old[i] = interp(R_old[i-1], R_old[i+1]); // For E_grav
            // U_old[i] = interp(U_old[i-1], U_old[i+1]);
        }
        else
        {
            // EDGES
            // P_old[i] = interp(P_old[i-1], P_old[i+1]);
            // Q_old[i] = interp(Q_old[i-1], Q_old[i+1]); // not necessary yet
            // Mass [i] = interp(Mass[i-1], Mass[i+1]);
        }
    }

    // enforce_boundary_conditions_Neumann(zones, U_old, P_old, Q_old);

    // ONLY DO THIS TO START
    U_old[0]        = 0;
    U_old[1]        = 0;
    U_old[zones-1]  = 0;
    U_old[zones-2]  = 0;


    U_new[0]        = 0;
    U_new[1]        = 0;
    U_new[zones-1]  = 0;
    U_new[zones-2]  = 0;

    R_new[0]        = R_old[0];
    R_new[1]        = R_old[1];
    R_new[zones-1]  = R_old[zones-1];
    R_new[zones-2]  = R_old[zones-2];

    V_new[0]        = V_old[0];
    V_new[zones-1]  = V_old[zones-1];



    //  Calculate initial energies
    E_kin  = calc_E_kin( Mass, U_old, zones);
    E_int  = calc_E_int( Mass, T_old, zones);
    E_grav = 0;
    if (with_gravity == true)
    {
        E_grav = calc_E_grav(Mass, M_int, R_old, zones);
    }

    E_tot  = E_kin + E_int + E_grav;

    Momentum = calc_Momentum(Mass, U_old, zones);






    // ========== PRINT INITIAL CONDITIONS ============ //
    

    char  info_fmt[] = "%10d %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le \n";
    char velos_fmt[] = "%10d %10d %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le \n";

    if (strcmp(IC, "restart")==0)
    {
         info_file_pointer       = fopen( info_filename, "a");
        velos_file_pointer       = fopen(velos_filename, "a");
    }
    else
    {
         info_file_pointer       = fopen( info_filename, "w");
        velos_file_pointer       = fopen(velos_filename, "w");
        fprintf(info_file_pointer, "%s %s %s %s %s %s %s %s %s %s",
            "# \t k  \t ",
            "  E_tot   \t\t ",
            "  E_grav  \t ",
            "  E_kin   \t\t ",
            "  E_int   \t ",
            "  Momentum   \t ",
            "  M_tot   \t\t ",
            "  delta_t \t ",
            "  time_current ",
            "\n");
        fprintf(velos_file_pointer, "%s %s %s %s %s %s %s %s %s %s %s %s",
            "# \t k  \t\t ",
            "  i  \t\t ",
            "  R_old [cm]  \t\t ",
            "  U_old [cm s^-1] \t ",
            "  1/V_old [g cm^-3]  \t ",
            "  T_old [K] \t\t ",
            "  Mass [g] \t\t",
            "  C_ad [cm s^-1]  \t ",
            "  E_old [ergs g^-1] \t ",
            "  P_old [dyne cm^-2] \t",
            "  Q_old [dyne cm^-2]",
            "\n");

        // PRINT ENERGY CONSERVATION INFORMATION
        fprintf(info_file_pointer, info_fmt, k_start, 
            E_tot, E_grav, E_kin, E_int, Momentum, 
            M_int[zones-2], 0., time_current);

        // PRINT VARIABLE INFORMATION
        for(i=0;i<zones;i++)
        {
            fprintf(velos_file_pointer, velos_fmt, 
                k_start, i, R_old[i], U_old[i], 1/V_old[i], T_old[i], Mass[i],
                C_ad[i], E_old[i], P_old[i], Q_old[i]);
        }
    }


    // ========== BEGIN RUN ============ //

    // double R_shock = 1.115 * pow(E_SN * M_SN * V_0,.2) * pow(time_current, .4);
    // START TIMESTEPS
    for (k=1; k<timesteps; ++k)
    {

        // R_shock = 1.115 * pow(E_SN * M_SN * V_0,.2) * pow(time_current, .4);
        // check if grid needs to be expanded -- improve this later
        // if (.8* R_old[zones-2] <  R_shock)
        // {
        //     printf("Need to expand zones! \n");
        //     printf("R_shock = %le \n", R_shock);
        //     printf("R_outer = %le \n", R_old[zones-2]);
        //     printf("k = %d \n", k);
        //     // if theoretical shock has made it 80% across zones, expand grid
        //     expand_grid(&zones, 
        //         &U_old, &U_new,
        //         &R_old, &R_new,
        //         &V_old, &V_new,
        //         &T_old, &T_new,
        //         &E_old, &E_new, &dE,
        //         &P_old, &P_new, 
        //         &Q_old, &Q_new,
        //         &C_ad,
        //         &Mass,  &M_int);

        //     printf("new 'zones' = %d \n", zones);
        //     printf("---------- \n");
        // }



        eta_var = 1; // If needed, will reduce eta_var to prevent negative energies
                   
        eta_var_continue = true;   // set eta_var_continue = false when iteration finished
        j = 0;     // Count how many times you've iterated eta_var

        while((eta_var_continue == true) && j<100)
        {
            i = 1;  // FIRST EDGE
            time_delta       = eta * eta_var * (R_old[i+2] - R_old[i]) 
                / sqrt( pow(C_ad[i+1],2) 
                        + ( (pow(U_old[i+2],2) + pow(U_old[i],2)) / 2) );

            if(!isfinite(time_delta))
            {
                printf("Error at innermost dt check \n");
                printf("k = %d \n", k);
                printf("(R_old[%d] - R_old[%d]) = %le \n", i+2,i,(R_old[i+2] - R_old[i]) );
                printf("R_old[%d] = %le \n", i+2, R_old[i+2]);
                printf("R_old[%d] = %le \n", i, R_old[i]);
                printf("denom = %le \n", sqrt( pow(C_ad[i+1],2) + ((pow(U_old[i+2],2) + pow(U_old[i],2)) / 2)));
                printf("pow(C_ad[i+1],2) = %le \n", pow(C_ad[i+1],2));
                printf("((pow(U_old[i+2],2) + pow(U_old[i],2)) / 2) = %le \n", ((pow(U_old[i+2],2) + pow(U_old[i],2)) / 2));
                return 0;
            }


            for (i=1; i<zones-3; i+=2)
            {
                // CHECK CFL CONDITION ACROSS EACH CELL
                time_delta_tmp   = eta * eta_var * (R_old[i+2] - R_old[i]) 
                    / sqrt( pow(C_ad[i+1],2) + ((pow(U_old[i+2],2) + pow(U_old[i],2)) / 2)); 


                if (time_delta_tmp < time_delta)
                {
                    time_delta = time_delta_tmp;
                    if(time_delta <= 0 || !isfinite(time_delta))
                    {
                        // to do: throw a proper error
                        printf("\n ==== ERROR ==== \n");
                        printf("delta t < 0 at i=%d\n", i);
                        printf("delta t = %le\n", time_delta_tmp);
                        printf("k = %d \n", k);
                        printf("eta = %lf \n", eta);
                        printf("eta_var = %le \n", eta_var);
                        printf("(R_old[i+2] - R_old[i]) = %le \n", (R_old[i+2] - R_old[i]) );
                        printf("R_old[i+2] = %le \n", R_old[i+2]);
                        printf("R_old[i] = %le \n", R_old[i]);
                        printf("denom = %le \n", sqrt( pow(C_ad[i+1],2) + ((pow(U_old[i+2],2) + pow(U_old[i],2)) / 2)));
                        printf("pow(C_ad[i+1],2) = %le \n", pow(C_ad[i+1],2));
                        printf("((pow(U_old[i+2],2) + pow(U_old[i],2)) / 2) = %le \n", ((pow(U_old[i+2],2) + pow(U_old[i],2)) / 2));
                        return 0;
                    }
                }
            }

            // Calculate du / dt at each EDGE
            for (i=3; i<zones-3; i+=2)
            {
                // NOT CONSERVATIVE!
                accel_thermo = 4 * M_PI * pow(R_old[i],2)
                    * ( (P_old[i-1]+Q_old[i-1]) - (P_old[i+1]+Q_old[i+1]) )
                    / Mass[i];
                accel_grav = 0;
                if (with_gravity == true)
                {
                    accel_grav = - G * M_int[i] / pow(R_old[1],2);
                }
                if (!isfinite(accel_grav))
                {
                    printf("ERROR in accel_grav = %le \n", accel_grav);
                    printf("k = %d \n", k);
                    printf("i = %d \n", i);
                    return 0;
                }

                if (!isfinite(accel_thermo))
                {
                    printf("ERROR in accel_thermo = %le \n", accel_thermo);
                    printf("k = %d \n", k);
                    printf("i = %d \n", i);
                    printf("P_old[%d-1] = %le \n", i, P_old[i-1]);
                    printf("P_old[%d+1] = %le \n", i, P_old[i+1]);
                    printf("Q_old[%d-1] = %le \n", i, Q_old[i-1]);
                    printf("Q_old[%d+1] = %le \n", i, Q_old[i+1]);

                    printf("Mass[%d] = %le \n", i-1, Mass[i-1]);
                    printf("Mass[%d] = %le \n", i, Mass[i]);
                    printf("Mass[%d] = %le \n", i+1, Mass[i+1]);
                    printf("Mass[%d] = %le \n", i+2, Mass[i+2]);
                    return 0;
                }

                U_new[i] = U_old[i] + time_delta * (accel_grav + accel_thermo);
            }

            // Calculate dr /dt at each EDGE
            for (i=3; i<zones-3; i+=2)
            {
                R_new[i] = R_old[i] + U_old[i] * time_delta;
            }


            // Calculate new densities in each ZONE
            for (i=2; i<zones-2; i+=2)
            {
                V_new[i] = (4./3) * M_PI
                    * ( pow(R_new[i+1],3) - pow(R_new[i-1],3) )
                    / Mass[i];
            }

            // Calculate E in each ZONE -- (internal energy)
            //      Using adiabatic condition, to integrate PdV heating within a zone
            //       - only first order accurate if cooling is included
            for (i=2; i<zones-2; i+=2)
            {
                // dE[i] = - (P_old[i]+Q_old[i]) * (V_new[i]-V_old[i]) ; // Approx
                dE[i] = - (P_old[i]+Q_old[i]) * pow(V_old[i], gamma)
                    * ( pow(V_new[i], 1-gamma) - pow(V_old[i], 1-gamma) ) / (1-gamma) ; // Exact
            }            

            if (with_cooling == true)
            {
                calc_cooling(V_old, T_old, U_old, dE, zones, metallicity, time_delta, my_units);
            }

            for (i=2; i<zones-2; i+=2)
            {
                E_new[i] = E_old[i] + dE[i];
            }


            // Calculate artificial viscosities in ZONES
            for (i=2; i<zones-2; i+=2)
            {
                if (V_old[i] > V_new[i])
                {
                    // Only add viscosity for compression

                    // Issues with interpolating this back to the zone edges?
                    Q_new[i] = (1/V_new[i]) * pow(a_visc,2)
                        * pow(U_new[i+1] - U_new[i-1], 2);    
                }
                else
                {
                    Q_new[i] = 0;
                }
            }


            eta_var_continue = false; // Change if we find a bad value
            for (i=2; i<zones-2; ++i)
            {
                if(i%2==0)
                {
                    // ZONES
                    if(E_new[i] < 0)
                    {
                        eta_var_continue = true;
                        // printf("using energy CFL at zone i=%d, time k=%d, iteration j=%d \n", i, k, j);
                    }
                    if((Q_new[i] < 0) || !isfinite(Q_new[i]))
                    {
                        eta_var_continue = true;
                        // printf("using Q CFL at zone i=%d, time k=%d, iteration j=%d \n", i, k, j);
                    }
                }
                else
                {
                    // EDGES
                    if (R_new[i] > R_new[i+2])
                    {
                        printf("ZONE CROSSING at i=%d and %d, j=%d \n", i, i+2, j);
                        eta_var_continue = true;
                        return 0;
                    }                    
                }
            }


            eta_var *= .1;
            ++j;
        }

        // printf("at k=%d, time_delta: %15.14le, time_current=%15.14le, iterations needed=%d \n", k, time_delta, time_current, j-1);
        if(!isfinite(time_delta))
        {
            printf("ERROR in time_delta! \n");
            return 0;
        }
        time_current += time_delta;

        // ========== UPDATE IMPLICITLY TRACKED VARIABLES ============ //


        // Calculate Pressure using Equation of State
        for (i=2; i<zones-2; i+=2)
        {
            P_new[i] = (gamma - 1) * E_new[i] / V_new[i];
        }

        // Calculate Temperature
        for (i=2; i<zones-2; i+=2)
        {
            T_new[i] = E_new[i] / c_V;
        }

        // Calculate sound speed in ZONES (for use with CFL check)
        for (i=2; i<zones-2; i+=2)
        {
            C_ad[i] = sqrt(gamma * P_new[i] * V_new[i]);
        }



        // ENFORCE BOUNDARY CONDITIONS
        // enforce_boundary_conditions(zones, 
        //     U_old, U_new,
        //     R_old, R_new,
        //     V_old, V_new,
        //     T_old, T_new,
        //     E_old, E_new,
        //     P_old, P_new,
        //     Q_old, Q_new );


        // // MOVE NEW VALUES -> OLD ARRAY
        swap(&U_old, &U_new);
        swap(&R_old, &R_new);
        swap(&V_old, &V_new);
        swap(&T_old, &T_new);
        swap(&E_old, &E_new);
        swap(&P_old, &P_new);
        swap(&Q_old, &Q_new);


        // interpolate values as necessary
        for(i=1; i<zones-1; ++i)
        {
            if(i%2==0)
            {
                // ZONES
                R_old[i] = interp(R_old[i-1], R_old[i+1]); // For E_grav
                // U_old[i] = interp(U_old[i-1], U_old[i+1]);
            }
            else
            {
                // EDGES
                // P_old[i] = interp(P_old[i-1], P_old[i+1]);
                // Q_old[i] = interp(Q_old[i-1], Q_old[i+1]); // not necessary yet
                // Mass [i] = interp(Mass[i-1], Mass[i+1]);
            }
        }



        // ========== PRINT VALUES ============ //


        if (k%n_print==0)
        {
            //  PRINT CURRENT VALUES
/*            printf("at k=%d, time_delta = %le, time_current=%le \n",
                k, time_delta, time_current);*/

            //  CALCULATE CURRENT ENERGIES
            E_kin  = calc_E_kin( Mass, U_old, zones);
            E_int  = calc_E_int( Mass, T_old, zones);
            E_grav = 0;
            if (with_gravity == true)
            {
                E_grav = calc_E_grav(Mass, M_int, R_old, zones);
            }
            E_tot  = E_kin + E_int + E_grav;

            Momentum = calc_Momentum(Mass, U_old, zones);

            // PRINT ENERGY CONSERVATION INFORMATION
            fprintf(info_file_pointer, info_fmt, k+k_start,
                E_tot, E_grav, E_kin, E_int, Momentum,
                M_int[zones-2], time_delta, time_current);

            // PRINT VARIABLE INFORMATION
            for(i=0;i<zones;++i)
            {
                fprintf(velos_file_pointer, velos_fmt, 
                    k+k_start, i, R_old[i], U_old[i], 1/V_old[i], T_old[i], Mass[i],
                    C_ad[i], E_old[i], P_old[i], Q_old[i]);
            }   
        }
     


        // ANYTHING ELSE?


    }

    // ========== CREATE A RESTART FILE ============ //

    FILE *mod_file_pointer; // for restarting runs
    mod_file_pointer    = fopen(mod_filename, "w");

    printf("Writing modfile to: %s \n", mod_filename);


    fprintf(mod_file_pointer, "%s \n", "# TIME (CURRENT) [s]");
    fprintf(mod_file_pointer, "%20.10le \n", time_current);
    fprintf(mod_file_pointer, "%s \n", " # TIMESTEP (k)");
    fprintf(mod_file_pointer, "%d \n", k+k_start);
    fprintf(mod_file_pointer, "%s \n", " # ZONES");
    fprintf(mod_file_pointer, "%d \n", zones);


    char mod_fmt[] = "%20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le \n";
    fprintf(mod_file_pointer, "%s %s %s %s %s %s %s %s %s",
        "# R_old [cm]  \t\t ",
        "  U_old [cm s^-1] \t ",
        "  V_old [g cm^-3]  \t ",
        "  T_old [K] \t ",
        "  E_old [ergs g^-1] \t ",
        "  P_old [dyne cm^-2] \t",
        "  Q_old [dyne cm^-2]",
        "  C_ad [cm s^-1]  \t ",
        "\n");

    for (i=0; i<zones; ++i)
    {
        fprintf(mod_file_pointer, mod_fmt, 
            R_old[i], U_old[i], V_old[i], T_old[i],
            E_old[i], P_old[i], Q_old[i], C_ad[i]); 
    }

    fclose(velos_file_pointer);
    fclose(info_file_pointer);
    fclose(mod_file_pointer);

    return 1;
};


double calc_E_kin(double Mass[], double U[], int zones )
{
    double E_kin = 0;
    int i;

    for(i=1; i<zones-1; i+=2)
    {
        E_kin += (1./2) * Mass[i] * pow(U[i], 2);
    }

    return E_kin;
};

double calc_E_int(double Mass[], double T[], int zones )
{
    double E_int = 0;
    int i;

    for(i=2; i<zones-2; i+=2)
    {
        E_int += Mass[i] * c_V * T[i];
    }

    return E_int;
};

double calc_E_grav(double Mass[], double M_int[], double R[], int zones)
{
    double E_grav = 0;
    int i;

    if (R[1]!=0)
    {
        E_grav = - (3./5) * G * pow(Mass[0], 2) / R[1];
    }

    for(i=2; i<zones-2; i+=2)
    {
        E_grav += - G * M_int[i-1] * Mass[i] / R[i];
    }

    return E_grav;
};

double calc_Momentum(double Mass[], double U[], int zones)
{
    double Momentum = 0;
    int i;

    for(i=1; i<zones-1; i+=2)
    {
        Momentum +=  Mass[i] * U[i];
    }

    return Momentum;    
}




