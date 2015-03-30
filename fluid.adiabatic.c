#include <stdio.h>
#include <stdbool.h>
#include <math.h>

double interp(double f_inner, double f_outer);

double calc_E_kin( double Mass[], double U[], int zones );
double calc_E_int( double Mass[], double T[], int zones, double c_V );
double calc_E_grav(double Mass[], double M_int[], double R[], int zones, double G );

int main()
{
    const int    timesteps   = 2000000; 
    const int    n_print     = timesteps / 1000; // print every n_steps
    // const int    timesteps   = 100000; 
    // const int    n_print     = 100; // print every n_steps

 
    double       E_int, E_kin, E_grav, E_tot;

    // all zone number variables must be ODD
    const int    zones       = 369; //remember to add 4 guard cells
    const int    zones_inner = 300;


    const double gamma       = 5./3;      // adiabatic index
    const double mu          = .67;      // mean molecular weight

    // const double pi          = 3.14159;
    const double G           = 6.67e-8;
    const double k_boltzmann = 1.38046e-16;
    const double m_proton    = 1.6598e-24;
    const double comp        = k_boltzmann / (m_proton * mu);

    const double R_gas       = 8.317e7;  // gas constant
    const double c_V         = (1./(gamma - 1)) * R_gas / mu; // specific heat capacity, const. volume


    const double AU          = 1.496e+13;
    const double pc          = 3.086e+18;
    const double R_sun       = 3.0e+11;
    const double M_sun       = 2.0e+33;


    // set these more properly later
    const double M_SN        = 3 * M_sun;
    const double E_SN        = 1e50 / M_SN; // Energy per unit mass
    

    const double R_total     = 1e2 * pc;  // total size, initially
    const double V_0         = 1. / (m_proton); // initial background density
    const double T_0         = 1e4;             // initial background temperature

    const double V_core      = 1. / (100 * m_proton);
    const double R_core      = 0; // radius of innermost zone boundary (unresolved) 
    const double M_core      = (4./3) * M_PI * pow(R_core,3) * V_core; // mass internal to R_core
    const double P_core      = 0;
    const double Q_core      = 0;
    const double E_core      = 0;
    const double T_core      = 0;
    const double C_ad_core   = sqrt(gamma * P_core * V_core);

    // const double V_background      = V_0;
    // const double T_background      = T_0;

    const double R_inner     = 1.5 * pc;   // size of beginning of sedov phase
    const double V_inner     = (4./3) * M_PI * (pow(R_inner,3) - pow(R_core,3)) / M_SN;
    const double T_inner     = E_SN / c_V;
    const double P_inner     = comp * T_inner / V_inner;



    const double a_visc           = 3;     // sets artificial viscosity strength
    const double eta              = 0.2;   // CFL number
          double eta_var          = 1.;    // Variable CFL knob
          bool   eta_var_continue = true;  // keep iterating to lower eta_var's?


    double       U_old[zones], U_new[zones];
    double       R_old[zones], R_new[zones];
    double       V_old[zones], V_new[zones];
    double       T_old[zones], T_new[zones];
    double       E_old[zones], E_new[zones];
    double       P_old[zones], P_new[zones];
    double       Q_old[zones], Q_new[zones];
    double       C_ad[ zones];
    double       Mass[ zones];
    double       M_int[zones];


    int          i,j,k;  // loop variables

    double       time_current = 0;

    // declare various temporary parameters
    double       time_delta_tmp, time_delta;
    double       accel_thermo, accel_grav;

    const bool   interactive = false;
    bool         from_previous, with_gravity;


    FILE   *in_file_pointer; // for restarting runs
    FILE *info_file_pointer, *velos_file_pointer, *out_file_pointer; // outputs
    char     in_filename[]   = "input";
    char   info_filename[]   = "info";
    char  velos_filename[]   = "velos";
    char    out_filename[]   = "output";
       in_file_pointer       = fopen(   in_filename, "r"); 
     info_file_pointer       = fopen( info_filename, "w");
    velos_file_pointer       = fopen(velos_filename, "w");
      out_file_pointer       = fopen(  out_filename, "w");
    char  info_fmt[] = "%10d %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le \n";
    char velos_fmt[] = "%10d %10d %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le %20.10le \n";
    fprintf(info_file_pointer, "%s %s %s %s %s %s %s %s %s",
        "# \t\t k  \t\t ",
        "  E_tot   \t\t ",
        "  E_grav  \t\t ",
        "  E_kin   \t\t ",
        "  E_int   \t\t ",
        "  M_tot   \t\t ",
        "  delta_t \t\t ",
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

    if (interactive == true)
    {
        int tmp; // for converting to boolean
        printf("Restart from model (1) or new simulation (0)?\n");
        scanf("%d", &tmp);
        from_previous = tmp;

        printf("Exclude gravity?  Yes (1) or no (0)\n");
        scanf("%d", &tmp);
        with_gravity = tmp;
    }
    else
    {
        from_previous = true;
        with_gravity  = false;
    }

    if (from_previous == true)
    {
        // READ DATA TO CONTINUE PREVIOUS RUN
        char buffer[1024];  // buffer for reading file
        printf("\nReading in from file: %s \n", in_filename);

        fgets(buffer, sizeof(buffer), in_file_pointer); // skip headers
        fscanf(in_file_pointer,"%lf \n", &time_current);
        // fgets(buffer, sizeof(buffer), in_file_pointer); // skip headers
        fgets(buffer, sizeof(buffer), in_file_pointer); // skip headers

        // printf("starting time: %lf [s]\n", time_current);
        for (i=0; i<zones; ++i)
        {
            fscanf(in_file_pointer, 
                "%lf %lf %lf %lf %lf %lf %lf %lf \n", 
                &(U_old[i]),
                &(R_old[i]),
                &V_old[i],
                &T_old[i],
                &E_old[i],
                &P_old[i],
                &Q_old[i],
                &C_ad [i]); 

            // printf("U_old[%d] = %lf \n", i, U_old[i]);
        }
        fclose(in_file_pointer);
    }
    else
    {
        // INITIALIZE DATA FOR A NEW RUN
        for (i=0; i<zones; i++)
        {
            // set defaults
            U_old[i] = 0;  // this might be overwritten for inner zones
            Q_old[i] = 0;

            if (i==0)
            {
                //  UNRESOLVED INNER CORE (zone)
                //      effectively acts as a guard cell
                Mass[ i] = M_core;

                V_old[i] = V_core;
                T_old[i] = T_core;
                E_old[i] = E_core; 
                P_old[i] = 0;           // get from neumann boundary conditions
                Q_old[i] = Q_core;
                C_ad[ i] = C_ad_core;
            }
            else if (i==1)
            {
                //  UNRESOLVED INNER CORE (boundary)
                //      effectively acts as a guard cell
                R_old[i] = R_core;

                // these are a little poorly defined for boundaries
                //      but it sets overall boundary conditions
                V_old[i] = V_core;
                T_old[i] = T_core;
                E_old[i] = E_core; 
                P_old[i] = 0;           // get from neumann boundary conditions
                Q_old[i] = Q_core;
                C_ad[ i] = C_ad_core;
            }
            else
            { 
            // Not inner guard cells:
                if(i<zones_inner)
                {
                    // inner pre-sedov zones
                    if (i%2 == 0)
                    {
                        //  ZONES
                        T_old[i] = T_inner;
   
                    }
                    else
                    {
                        //  BOUNDARIES
                        R_old[i] = i / ((double) zones_inner) * R_inner;
   
                    }
                    // variables which don't particularly care 
                    //  if they're in the zone or boundary
                    V_old[i] = V_inner;
                    T_old[i] = T_inner;
                    E_old[i] = E_SN;
                    P_old[i] = P_inner;
                    C_ad[ i] = sqrt(gamma * P_old[i] * V_old[i]);
                }
                else
                {
                    // uniform background
                    if(i%2 == 0)
                    {
                        // background ZONES
                            // information set below.
                            //      V, T, etc are the same between zones
                    }
                    else
                    {
                        // background BOUNDARIES
                        R_old[i] = R_inner
                            + ( (i-zones_inner) / ((double) zones-zones_inner)
                                * (R_total-R_inner) );
                    }

                    // background variables which could be zones or boundaries
                    V_old[i] = V_0;
                    T_old[i] = T_0;
                    E_old[i] = c_V * T_0;
                    P_old[i] = comp * T_0 / V_0;
                    C_ad[ i] = sqrt(gamma * P_old[i] * V_old[i]);
                }                
            }
        }
    }

    //  Calculate Mass of each zone
    for(i=0; i<zones-2; i+=2)
    {
        if (i==0)
        {
            Mass[i] = (4./3) * M_PI * pow(R_old[i+1],3) / V_old[i];
        }
        else
        {
            Mass[i] = (4./3) * M_PI * (pow(R_old[i+1],3) - pow(R_old[i-1],3))
                / V_old[i];
        }
    }

    //  Calculate Mass internal to each boundary
    M_int[1] = Mass[0];
    for(i=3; i<zones-1; i+=2)
    {
        M_int[i] = M_int[i-2] + Mass[i-1];
    }


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
            // BOUNDARIES
            // P_old[i] = interp(P_old[i-1], P_old[i+1]);
            // Q_old[i] = interp(Q_old[i-1], Q_old[i+1]); // not necessary yet
            Mass [i] = interp(Mass[i-1], Mass[i+1]);
        }
    }

        // UPDATE GUARD CELLS + BOUNDARY CONDITIONS
        //      INSIDE
        //          Dirichlet: 
        //              R = R_core, U = 0
        //              P = 0  // So provides no momentum / energy. E.g. hard wall
        //              Shouldn't need to be updated,
        //              if I did things right.
        //              (Subject to change depending on sub-grid model)
        //          Neumann:
        //              Maybe pressure, to avoid strange momentum fluxes
        P_old[0] = P_old[2];
        P_old[1] = P_old[2];
        Q_old[0] = Q_old[2];
        Q_old[1] = Q_old[2];
        U_old[0] = 0; // should be implicitly enforced by neumann on P,Q 
        U_old[1] = 0; 


        //      OUTSIDE
        //          Dirichlet:
        //              Temperature, Density, u
        //          Neumann:
        //              P - dP / dr = 0 at i=zones-1
        //                  - prevents momentum flux until the shock gets there
        //                    but by that point the simulation breaks anyway
        //                  - so this is basically dirichlet?
        //                      - but Neumann ensures no momentum flux,
        //                        even if the initial condition was poorly chosen

        P_old[zones-1] = P_old[zones-3] 
            + (P_old[zones-3]-P_old[zones-5])*(Mass[i]/Mass[i-2]);
        // P_old[zones-2] = P_old[zones-3];
        Q_old[zones-1] = Q_old[zones-3] 
            + (Q_old[zones-3]-Q_old[zones-5])*(Mass[i]/Mass[i-2]);
        // Q_old[zones-2] = Q_old[zones-3];
        U_old[zones-1] = 0;
        U_old[zones-1] = 0;



    //  Calculate initial energies
    E_kin  = calc_E_kin( Mass, U_old, zones);
    E_int  = calc_E_int( Mass, T_old, zones, c_V);
    E_grav = 0;
    if (with_gravity == true)
    {
        E_grav = calc_E_grav(Mass, M_int, R_old, zones, G);
    }

    E_tot  = E_kin + E_int + E_grav;

    // PRINT ENERGY CONSERVATION INFORMATION
    fprintf(info_file_pointer, info_fmt, 0, E_tot, E_grav, E_kin, E_int, M_int[zones-2], 0., time_current);

    // PRINT VARIABLE INFORMATION
    for(i=0;i<zones;i++)
    {
        fprintf(velos_file_pointer, velos_fmt, 
            0, i, R_old[i], U_old[i], 1/V_old[i], T_old[i], Mass[i],
            C_ad[i], E_old[i], P_old[i], Q_old[i]);
    }



    // START TIMESTEPS
    for (k=1; k<timesteps; ++k)
    {

        eta_var = 1; // If needed, will reduce eta_var to prevent negative energies
                   
        eta_var_continue = true;   // set eta_var_continue = false when iteration finished
        j = 0;     // Count how many times you've iterated eta_var

        while((eta_var_continue == true) && j<100)
        {
            i = 1;  // FIRST BOUNDARY
            time_delta       = eta * eta_var * (R_old[i+2] - R_old[i]) 
                / sqrt( pow(C_ad[i+1],2) 
                        + ( (pow(U_old[i+2],2) + pow(U_old[i],2)) / 2) );

            if(!isfinite(time_delta))
            {
                printf("Error at innermost dt check \n");
            }


            for (i=3; i<zones-3; i+=2)
            {
                // CHECK CFL CONDITION ACROSS EACH CELL
                time_delta_tmp   = eta * eta_var * (R_old[i+2] - R_old[i]) 
                    / sqrt( pow(C_ad[i+1],2) + (pow(U_old[i+2],2) + pow(U_old[i],2)) / 2); 

                // time_delta_tmp_E = 1 * eta * E_old[i+1] * (R_old[i+2] - R_old[i])
                //     / abs( (P_old[i+2]*U_old[i+2]) - (P_old[i]*U_old[i]) );

                // if (time_delta_tmp_E < time_delta)
                // {
                //     printf("\n Using an energy CFL at k=%d, i=%d", k, i);
                //     time_delta = time_delta_tmp_E;
                // }

                if (time_delta_tmp < time_delta)
                {
                    time_delta = time_delta_tmp;
                    if(time_delta <= 0 || !isfinite(time_delta))
                    {
                        // to do: throw a proper error
                        printf("\n ==== ERROR ==== \n");
                        printf("delta t < 0 at i=%d\n", i);
                        printf("delta t = %le\n", time_delta_tmp);
                    }
                }
            }

            // Calculate du / dt at each BOUNDARY
            for (i=1; i<zones-1; i+=2)
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

                // if(i==1)
                // {
                //     printf("accel_thermo: %le \n", accel_thermo);
                //     printf("Mass[i]: %le \n", Mass[i]);
                // }

                U_new[i] = U_old[i] + time_delta * (accel_grav + accel_thermo);
            }

            // Calculate dr /dt at each BOUNDARY
            for (i=1; i<zones-1; i+=2)
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
                E_new[i] = E_old[i]
                    - (P_old[i]+Q_old[i]) * (V_new[i]-V_old[i]) ;
            }

            // Calculate artificial viscosities in ZONES
            for (i=2; i<zones-2; i+=2)
            {
                if (V_old[i] > V_new[i])
                {
                    // Only add viscosity for compression

                    // Issues with interpolating this back to the zone boundaries?
                    Q_new[i] = (1/V_new[i]) * pow(a_visc,2)
                        * pow(U_new[i+1] - U_new[i-1], 2);    
                    // printf("adding Q = %le to zone %d at timestep %d \n",Q_new[i], i, k);
                    // if(i==2)
                    // {
                    //     printf("1/V_new[i] = %le \n", (1/V_new[i]));
                    //     printf("pow(a_visc,2) = %le \n", pow(a_visc,2));
                    //     printf("U_new[i+1] - U_new[i-1] = %le \n", U_new[i+1] - U_new[i-1]);
                    //     printf("pow(U_new[i+1] - U_new[i-1], 2) = %le \n", pow(U_new[i+1] - U_new[i-1], 2));
                    //     printf("U_new[i+1] %le \n", U_new[i+1]);
                    //     printf("U_new[i-1] %le \n", U_new[i-1]);
                    // }
                }
                else
                {
                    Q_new[i] = 0;
                }
            }


            eta_var_continue = false; // Change if we find a bad value
            for (i=1; i<zones-2; ++i)
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
                    // BOUNDARIES
                    if (R_new[i] > R_new[i+2])
                    {
                        printf("ZONE CROSSING at i=%d and %d, j=%d", i, i+2, j);
                        eta_var_continue = true;
                    }                    
                }
            }


            eta_var *= .1;
            ++j;
        }
        // printf("at k=%d, time_delta: %15.14le, time_current=%15.14le, iterations needed=%d \n", k, time_delta, time_current, j-1);
        if(!isfinite(time_delta))
        {
            printf("ERROR! \n");
        }
        time_current += time_delta;


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

        // Interpolate necessary values
        for (i=1; i<zones-1; ++i)
        {
            if(i%2==0)
            {
                // ZONES
                R_new[i] = (R_new[i-1] + R_new[i+1]) / 2; // for E_grav
                // U_new[i] = (U_new[i-1] + U_new[i+1]) / 2;
            }
            else
            {
                // BOUNDARIES
                V_new[i] = (V_new[i-1] + V_new[i+1]) / 2; // not necessary, but gets rid of nan's in velos
                // T_new[i] = (T_new[i-1] + T_new[i+1]) / 2;
                // E_new[i] = (E_new[i-1] + E_new[i+1]) / 2;
                // P_new[i] = (P_new[i-1] + P_new[i+1]) / 2;
                // Q_new[i] = (Q_new[i-1] + Q_new[i+1]) / 2;
                // C_ad [i] = (C_ad [i-1] + C_ad [i+1]) / 2;
            }
        }


        // MOVE NEW VALUES -> OLD ARRAY
        for (i=2; i<zones-2; ++i)
        {
            U_old[i] = U_new[i];
            R_old[i] = R_new[i];
            V_old[i] = V_new[i];
            T_old[i] = T_new[i];
            E_old[i] = E_new[i];
            P_old[i] = P_new[i];
            Q_old[i] = Q_new[i];
        }


        // UPDATE GUARD CELLS + BOUNDARY CONDITIONS
        //      INSIDE
        //          Dirichlet: 
        //              R = R_core, U = 0
        //              P = 0  // So provides no momentum / energy. E.g. hard wall
        //              Shouldn't need to be updated,
        //              if I did things right.
        //              (Subject to change depending on sub-grid model)
        //          Neumann:
        //              Maybe pressure, to avoid strange momentum fluxes
        P_old[0] = P_old[2];
        P_old[1] = P_old[2];
        Q_old[0] = Q_old[2];
        Q_old[1] = Q_old[2];
        U_old[0] = 0; // should be implicitly enforced by neumann on P,Q 
        U_old[1] = 0; 


        //      OUTSIDE
        //          Dirichlet:
        //              Temperature, Density, u
        //          Neumann:
        //              P - dP / dr = 0 at i=zones-1
        //                  - prevents momentum flux until the shock gets there
        //                    but by that point the simulation breaks anyway
        //                  - so this is basically dirichlet?
        //                      - but Neumann ensures no momentum flux,
        //                        even if the initial condition was poorly chosen
        //          Other(?):
        //              dP/dm - match last boundary derivative to previous boundary


        P_old[zones-1] = P_old[zones-3] 
            + (P_old[zones-3]-P_old[zones-5])*(Mass[i]/Mass[i-2]);
        // P_old[zones-2] = P_old[zones-3];
        Q_old[zones-1] = Q_old[zones-3] 
            + (Q_old[zones-3]-Q_old[zones-5])*(Mass[i]/Mass[i-2]);
        // Q_old[zones-2] = Q_old[zones-3];
        U_old[zones-1] = 0;
        U_old[zones-1] = 0;

        if (k%n_print==0)
        {
            //  PRINT CURRENT VALUES
/*            printf("at k=%d, time_delta = %le, time_current=%le \n",
                k, time_delta, time_current);*/

            //  CALCULATE CURRENT ENERGIES
            E_kin  = calc_E_kin( Mass, U_old, zones);
            E_int  = calc_E_int( Mass, T_old, zones, c_V);
            E_grav = 0;
            if (with_gravity == true)
            {
                E_grav = calc_E_grav(Mass, M_int, R_old, zones, G);
            }
            E_tot  = E_kin + E_int + E_grav;

            // PRINT ENERGY CONSERVATION INFORMATION
            fprintf(info_file_pointer, info_fmt, k, E_tot, E_grav, E_kin, E_int, M_int[zones-2], time_delta, time_current);

            // PRINT VARIABLE INFORMATION
            for(i=0;i<zones;++i)
            {
                fprintf(velos_file_pointer, velos_fmt, 
                    k, i, R_old[i], U_old[i], 1/V_old[i], T_old[i], Mass[i],
                    C_ad[i], E_old[i], P_old[i], Q_old[i]);
            }   
        }
     


        // ANYTHING ELSE?


    }


    return 0;
};

double interp(double f_inner, double f_outer)
{
    // For interpolating variables between zones and boundaries
    // 
    // General interpolates a value f_middle
    //      given f_inner, f_outer
    // 
    // This provides an easy way to consistently change the interpolation scheme

    double f_middle;
 
    f_middle = (f_inner + f_outer) / 2;

    return f_middle;

}

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

double calc_E_int(double Mass[], double T[], int zones, double c_V )
{
    double E_int = 0;
    int i;

    for(i=2; i<zones-2; i+=2)
    {
        E_int += Mass[i] * c_V * T[i];
    }

    return E_int;
};

double calc_E_grav(double Mass[], double M_int[], double R[], int zones, double G)
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
