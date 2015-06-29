#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

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
int    zones       = (2*1024) - 1; // including guard cells (remember to add 4 guard cells)

char   IC[80]      = "rt1d_comparison";

const bool   with_cooling = false;
const bool   with_gravity = false;

const bool   Q_tensor         = true;
const bool   Q_tensor_whalen  = true; // if Q_tensor == true, Q_tensor_whalen == false, then use Schulz
const bool   with_conduction  = true;

int main()
{

    if (zones%2 == 0)
    {
        zones += 1;
    }    

    int timesteps, n_print;

    timesteps   = 1000000; 
    // timesteps   = 4*10000; 
    // timesteps   = 100;
    // timesteps=2;

    n_print     = timesteps / 100; // print every n_steps
    if (n_print <= 0)
    {
        n_print = 1;
    }

    const double metallicity = grackle_data.SolarMetalFractionByMass;

          double a_visc_0         = 2;     // sets artificial viscosity strength
          double a_visc_1         = 1./5;     // sets artificial viscosity strength
          double h_0              = .005*4.;     // sets artificial conduction strength // See Noh 1987
          double h_1              = .005*1.;     // sets artificial conduction strength
    const double eta              = 0.2;   // CFL number
          double eta_var          = 1.;    // Variable CFL knob
          bool   eta_var_continue = true;  // keep iterating to lower eta_var's?

    if (strcmp(IC, "rt1d_comparison")==0)
    {
        h_0 = 1*4.;
        h_1 = 10*1.;
    }


    double       E_int, E_kin, E_grav, E_tot, Momentum; 

    double *U_old, *U_new;
    double *R_old, *R_new;
    double *V_old, *V_new;
    double *T_old, *T_new;
    double *E_old, *E_new, *dE;
    double *dE_cool;
    dE_cool = malloc(zones * sizeof(double));
    double *dE_conduction;
    dE_conduction = malloc(zones * sizeof(double));
    double *dE_adiabatic;
    dE_adiabatic = malloc(zones * sizeof(double));    
    double *P_old, *P_new;
    double *Q_old, *Q_new;
    double *H_old, *H_new;
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
    H_old = malloc(zones * sizeof(double));
    H_new = malloc(zones * sizeof(double));
    C_ad  = malloc(zones * sizeof(double));
    Mass  = malloc(zones * sizeof(double));
    M_int = malloc(zones * sizeof(double));

    int          i,j,k;  // loop variables
    int          k_start = 0; // start at timestep 0, unless continuing a run

    double       time_current = 0;

    // declare various temporary parameters
    double       time_delta_tmp, time_delta;


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
    printf("Using ICs: %s \n", IC);

    if (strcmp(IC, "sedov")==0)
    {
        time_current = sedov(zones, U_old, R_old, V_old, T_old, 
            E_old, P_old, Q_old, H_old, C_ad); 
    }
    else if (strcmp(IC, "restart")==0)
    {
        time_current = restart(&zones, &k_start, U_old, R_old, V_old, T_old, 
            E_old, P_old, Q_old, H_old, C_ad, mod_filename); 
        printf("k_start = %d \n", k_start);
    }
    else if (strcmp(IC, "new_blast")==0)
    {
        // INITIALIZE DATA FOR A NEW RUN
        // Only add energy to innermost zone
        time_current = new_blast(zones, U_old, R_old, V_old, T_old,
            E_old, P_old, Q_old, H_old, C_ad);
    }
    else if (strcmp(IC, "new_blast_spread")==0)
    {
        // INITIALIZE DATA FOR A NEW RUN
        // Only add energy to innermost zone
        time_current = new_blast_spread(zones, U_old, R_old, V_old, T_old,
            E_old, P_old, Q_old, H_old, C_ad);
    }
    else if (strcmp(IC, "rt1d_comparison")==0)
    {
        // INITIALIZE DATA FOR A NEW RUN
        // Only add energy to innermost zone
        time_current = rt1d_comparison(zones, U_old, R_old, V_old, T_old,
            E_old, P_old, Q_old, H_old, C_ad);
    }
    else
    {
        printf("ERROR: Please choose an initial condition.\n");        
        printf("Options: 'sedov', 'restart', or 'new_blast'.\n"); 
        return 0;       
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
            U_old[i] = interp(U_old[i-1], U_old[i+1]);
        }
        else
        {
            // EDGES
            // P_old[i] = interp(P_old[i-1], P_old[i+1]);
            Q_old[i] = interp(Q_old[i-1], Q_old[i+1]); // needed for conduction
            C_ad[i]  = interp(C_ad[i-1], C_ad[i+1]); // needed for conduction
            // Mass [i] = interp(Mass[i-1], Mass[i+1]);
        }
    }

    enforce_boundary_conditions_initial(zones, 
        U_old, U_new,
        R_old, R_new,
        V_old, V_new,
        T_old, T_new,
        E_old, E_new,
        P_old, P_new,
        Q_old, Q_new,
        H_old, H_new );


    //  Calculate initial global parameters
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
    for (i=0; i<zones; ++i)
    {
        dE[i] = 0;
    }

    char  info_fmt[] = "%10d %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le \n";
    char velos_fmt[] = "%10d %10d %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le \n";

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
        fprintf(velos_file_pointer, "%s %s %s %s %s %s %s %s %s %s %s %s %s",
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
            "  dE [ergs g^-1] \t ",
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
                C_ad[i], E_old[i], P_old[i], Q_old[i], dE[i]);
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
        //         &H_old, &H_new,
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
                double accel_thermo, accel_grav;

                // NOT CONSERVATIVE!
                accel_thermo = 4 * M_PI * pow(R_old[i],2)
                    * ( (P_old[i-1]+Q_old[i-1]) - (P_old[i+1]+Q_old[i+1]) )
                    / Mass[i]; 


                if (Q_tensor_whalen == true)
                {
                    // Noh 1987's Eq. 10.6
                    accel_thermo = accel_thermo - (3 * Q_old[i] * V_old[i] / R_old[i]); 
                }


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

            // Calculate dr /dt at each EDGE f
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
                
                // dE_adiabatic[i] = E_old[i] * (1 - pow(V_old[i]/V_new[i], gamma-1));
                // dE[i]           = dE_adiabatic[i];
                // dE[i]           = dE[i] + Q_old[i] * (V_old[i] - V_new[i]);


                // dE_adiabatic[i] = - (P_old[i]+Q_old[i]) * pow(V_old[i], gamma)
                    // * ( pow(V_new[i], 1-gamma) - pow(V_old[i], 1-gamma) ) / (1-gamma) ; // Exact

                dE_adiabatic[i] = - (P_old[i]+Q_old[i]) * (V_new[i]-V_old[i]) ; // Approx
                
                dE[i] = dE_adiabatic[i];



                if (Q_tensor_whalen == true)
                {
                    // DISABLED since it doesn't appear to conserve energy? Not sure what the problem is
                    // Noh 1987's Euqation 10.6 
                    dE[i] = dE[i] + ((3 * U_old[i] * Q_old[i] * V_old[i] / R_old[i]) * time_delta);
                }

                if (with_conduction == true)
                {
                    dE_conduction[i] = 3 * (pow(R_old[i+1],2)*H_old[i+1] - pow(R_old[i-1],2)*H_old[i-1]) * time_delta / Mass[i];

                    // if ((k==45775) && (i==10) )
                    if (false)
                    {
                        printf("dE adiabatic = %le \n", - (P_old[i]+Q_old[i]) * (V_new[i]-V_old[i]));
                        printf("dE viscous = %le \n", (3 * U_old[i] * Q_old[i] / R_old[i]));
                        printf("dE_conduction = %le \n", dE_conduction[i]);
                        printf("E_old[%d] = %le \n", i, E_old[i]);
                    }
                    dE[i] = dE[i] + 3 * (pow(R_old[i+1],2)*H_old[i+1] - pow(R_old[i-1],2)*H_old[i-1]) * time_delta / Mass[i]; 
                }
            }            

            if (with_cooling == true)
            {

                for (i=0; i<zones; ++i)
                {
                    dE_cool[i] = 0;
                }
                calc_cooling(V_old, E_old, U_old, dE_cool, zones, metallicity, time_delta, my_units);
                for(i=0; i<zones; ++i)
                {
                    dE[i] += dE_cool[i];
                }
            }

            for (i=2; i<zones-2; i+=2)
            {
                E_new[i] = E_old[i] + dE[i];
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
                        printf("using energy CFL at zone i=%d, time k=%d, iteration j=%d \n", i, k, j);
                        printf("E_old[%d] = %le \n", i, E_old[i]);
                        printf("T_old[%d] = %lf \n", i, T_old[i]);
                        printf("dE_cool = %le \n", dE_cool[i]);
                        printf("dE_conduction = %le \n", dE_conduction[i]);
                        printf("dE_adiabatic = %le \n", dE_adiabatic[i]);
                        printf("U_old[%d] = %le \n", i-1, U_old[i-1]);
                        printf("U_old[%d] = %le \n", i+1, U_old[i+1]);
                        printf("V_old[%d] = %le \n", i, V_old[i]);
                        printf("V_new[%d] = %le \n", i, V_new[i]);
                        return;
                    }
                    // if((Q_new[i] < 0) || !isfinite(Q_new[i]))
                    // {
                    //     eta_var_continue = true;
                    //     printf("using Q CFL at zone i=%d, time k=%d, iteration j=%d \n", i, k, j);
                    // }
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


        // Calculate artificial viscosities in ZONES
        for (i=2; i<zones-2; i+=2)
        {
            double delta_U, delta_R, U_avg, R_avg;
            delta_U = U_new[i+1] - U_new[i-1];
            if (Q_tensor == true)
            {
                if (Q_tensor_whalen == true)
                {
                    // Whalen tensor viscosity
                    if (delta_U < 0) // check using velocity gradient
                    {
                        delta_R = R_new[i+1] - R_new[i-1];

                        U_avg = interp(U_new[i+1], U_new[i-1]);
                        R_avg = interp(R_new[i+1], R_new[i-1]);

                        Q_new[i] = pow(a_visc_0,2) * (1/V_new[i])    // Whalen's tensor
                            * delta_U * (delta_U - (U_avg * delta_R / R_avg)) ;   
                    }
                    else
                    {
                        // not compressing
                        Q_new[i] = 0;
                    }
                }
                else
                {
                    // Schulz's artificial visocity tensory (Noh 1987, sec 9)
                    if (delta_U < 0)
                    {
                        // Noh 1987, eq 9.3
                        double U_rr; // second derivative 

                        Q_new[i] = 0;
                    }
                    else{
                        // not compressing
                        Q_new[i] = 0;
                    }
                }
            }
            else
            {
                // von Neumann viscosity
                if (V_old[i] > V_new[i]) // check compression using volumes
                {
                    Q_new[i] = pow(a_visc_0,2) * (1/V_new[i]) * pow(delta_U, 2)
                        - a_visc_1 * (1/V_new[i]) * C_ad[i] * delta_U;    
                }
                else
                {
                    // not compressing
                    Q_new[i] = 0;
                }
            }

            if (Q_new[i] < 0)
            {
                Q_new[i] = 0;
            }
        }


        // Calculate artificial conduction in ZONES
        for (i=3; i<zones-3; i+=2)
        {
            if (with_conduction == true)
            {
                if (((Q_new[i+1] > 1e-20)) || (Q_new[i-1] > 1e-20))
                {

                    double delta_U_L        = fabs(fmin(U_new[i  ] - U_new[i-2], 0));
                    double delta_U_R        = fabs(fmin(U_new[i+2] - U_new[i  ], 0));

                    double rho_delta_U_L    = (1/V_new[i-1]) * delta_U_L; 
                    double rho_delta_U_R    = (1/V_new[i+1]) * delta_U_R;

                    double rho_delta_U_avg  = harmonic_mean(rho_delta_U_L, rho_delta_U_R); 

                    if((rho_delta_U_L < 0) || (rho_delta_U_R < 0) || (rho_delta_U_avg < 0))
                    {
                        printf("Error in rho_delta_U's \n");
                        printf("k = %d, i = %i \n", k, i);
                        printf("rho_delta_U_L = %le \n", rho_delta_U_L);
                        printf("rho_delta_U_R = %le \n", rho_delta_U_R);
                        printf("rho_delta_U_avg = %le \n", rho_delta_U_avg); 
                    }

                    double rho_C_ad_L   = (1/V_new[i-1]) * C_ad[i-1]; 
                    double rho_C_ad_R   = (1/V_new[i+1]) * C_ad[i+1]; 

                    double rho_C_ad_avg = harmonic_mean(rho_C_ad_L, rho_C_ad_R);

                    double delta_E      =    E_new[i+1] - E_new[i-1]; // SPATIAL delta_E, not temporal dE


                    // See Noh (1987), Eqs. 2.3, 3.9
                    H_new[i] = pow(h_0,2) * rho_delta_U_avg * delta_E
                        + h_1 * rho_C_ad_avg * delta_E;

                    double term_0 = pow(h_0,2) * rho_delta_U_avg * delta_E;
                    double term_1 = h_1 * rho_C_ad_avg * delta_E;

                    if (!isfinite(H_new[i]))
                    {
                        printf("ERROR in H_new at i=%d, k=%d; H_new = %le \n", i,k, H_new[i]);
                        printf("Q_new[%d] = %le \n", i-1, Q_new[i-1]);
                        printf("Q_new[%d] = %le \n", i+1, Q_new[i+1]);
                        printf("h_0 = %le \n", h_0);
                        printf("(2*rho_delta_U_L * rho_delta_U_R) = %le \n", (2*rho_delta_U_L * rho_delta_U_R));
                        printf("(  rho_delta_U_L + rho_delta_U_R) = %le \n", (  rho_delta_U_L + rho_delta_U_R));
                        printf(" rho_delta_U_L = %le \n", rho_delta_U_L);
                        printf(" rho_delta_U_R = %le \n", rho_delta_U_R);
                        printf(" delta_E = %le \n", delta_E);
                        return;
                    }
                }
                else
                {
                    // not compressing on either side
                    H_new[i] = 0;
                }
            }
            else
            {
                // not using conduction
                H_new[i] = 0;
            }
        }

        // ENFORCE BOUNDARY CONDITIONS
        enforce_boundary_conditions(zones, 
            U_old, U_new,
            R_old, R_new,
            V_old, V_new,
            T_old, T_new,
            E_old, E_new,
            P_old, P_new,
            Q_old, Q_new,
            H_old, H_new );


        // // MOVE NEW VALUES -> OLD ARRAY
        swap(&U_old, &U_new);
        swap(&R_old, &R_new);
        swap(&V_old, &V_new);
        swap(&T_old, &T_new);
        swap(&E_old, &E_new);
        swap(&P_old, &P_new);
        swap(&Q_old, &Q_new);
        swap(&H_old, &H_new);


        // interpolate values as necessary
        for(i=1; i<zones-1; ++i)
        {
            if(i%2==0)
            {
                // ZONES
                R_old[i] = interp(R_old[i-1], R_old[i+1]); // For E_grav
                U_old[i] = interp(U_old[i-1], U_old[i+1]);
            }
            else
            {
                // EDGES
                // P_old[i] = interp(P_old[i-1], P_old[i+1]);
                // Q_old[i] = interp(Q_old[i-1], Q_old[i+1]); // needed for conduction
                C_ad[i]  = interp(C_ad[i-1], C_ad[i+1]); // needed for conduction
                // Mass [i] = interp(Mass[i-1], Mass[i+1]);
            }
        }



        // ========== PRINT VALUES ============ //


        if (k%n_print==0)
        {
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
                    C_ad[i], E_old[i], P_old[i], Q_old[i], dE_cool[i]);
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


    char mod_fmt[] = "%20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le \n";
    fprintf(mod_file_pointer, "%s %s %s %s %s %s %s %s %s %s",
        "# R_old [cm]  \t\t ",
        "  U_old [cm s^-1] \t ",
        "  V_old [g cm^-3]  \t ",
        "  T_old [K] \t ",
        "  E_old [ergs g^-1] \t ",
        "  P_old [dyne cm^-2] \t",
        "  Q_old [dyne cm^-2]",
        "  H_old [dyne cm^-2]",
        "  C_ad [cm s^-1]  \t ",
        "\n");

    for (i=0; i<zones; ++i)
    {
        fprintf(mod_file_pointer, mod_fmt, 
            R_old[i], U_old[i], V_old[i], T_old[i],
            E_old[i], P_old[i], Q_old[i], H_old[i], C_ad[i]); 
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

