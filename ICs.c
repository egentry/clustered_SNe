
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <constants.h> // defines physical constants

#include <ICs.h>


double sedov(int zones, double U[], double R[], double V[], double T[], 
            double E[], double P[], double Q[], double H[], double C_ad[])
{
    double time_current; 
    double R_blast;

    R_blast = .01 * pc;
    R_blast = pc;

    int sedov_zones = zones / 20;
    int outer_zones = zones - sedov_zones;
    char command[80];
    int ret;
    ret = sprintf(command, "python dimensionalize_sedov_initial.py %d %d %le %le", 
        sedov_zones, outer_zones, R_blast, R_total / 1 );
    printf("create sedov IC using command: %s \n", command);
    ret = system(command);


    FILE *in_file_pointer; // for restarting runs
    char  in_filename[]      = "input";
          in_file_pointer    = fopen(in_filename, "r"); 



    // READ DATA TO CONTINUE PREVIOUS RUN
    char buffer[1024];  // buffer for reading file
    printf("\nReading in from file: %s \n", in_filename);

    fgets(buffer, sizeof(buffer), in_file_pointer); // skip headers
    fscanf(in_file_pointer,"%lf \n", &time_current);
    fgets(buffer, sizeof(buffer), in_file_pointer); // skip headers

    printf("starting time: %lf [yr]\n", time_current/ yr);
    int i;
    for (i=0; i<zones; ++i)
    {
        fscanf(in_file_pointer, 
            "%lf %lf %lf %lf %lf %lf %lf %lf %lf \n", 
            &(U[i]),
            &(R[i]),
            &V[i],
            &T[i],
            &E[i],
            &P[i],
            &Q[i],
            &H[i],
            &C_ad[i]); 

    }
    fclose(in_file_pointer);

    return time_current;
}


double new_blast(int zones,
    double U[], double R[], double V[], double T[],
    double E[], double P[], double Q[], double H[], double C_ad[])
{

    printf("Initializing blast at innermost cell. \n");
    double time_current;
    int i;

    int inner_zones = zones / 10;
    inner_zones = 0;
    int outer_zones = zones - inner_zones;
    // double R_inner = R_total / 100;

    // INITIALIZE DATA FOR A NEW RUN
    for (i=0; i<zones; i++)
    {
        // set defaults
        U[i] = 0;
        
        if (i < inner_zones)
        {
            R[i] = i * R_inner / inner_zones;
        }
        else
        {
            R[i] = R_inner
                + (i-inner_zones) * (R_total - R_inner) / (zones - inner_zones);
        }

        // // POWERLAW SPACING
        // double power_law_index = 1.1;
        // double C = R_total / pow(zones-2, power_law_index);
        // R[i] = C * pow(i, power_law_index);

        V[i] = V_0;
        T[i] = T_0;
        E[i] = c_V * T[i];
        P[i] = (gamma - 1) * E[i] / V[i];
        Q[i] = 0;
        H[i] = 0;
        C_ad[i] = sqrt(gamma * P[i] * V[i]);
    }

    i = 2;
    // innermost real zone
    double M;
    V[i] = (4.*M_PI/3.) * (pow(R[i+1],3) - pow(R[i-1],3)) / M_SN;
    // V[i] = V_0;
    M = (1./V[i]) * (4.*M_PI/3.) * (pow(R[i+1],3) - pow(R[i-1],3));
    E[i] = E_SN * M_SN / M; // PER UNIT MASS!
    T[i] = E[i] / c_V;
    P[i] = (gamma - 1) * E[i] / V[i];
    C_ad[i] = sqrt(gamma * P[i] * V[i]);

    time_current = 0;
    return time_current;
}

double new_blast_spread(int zones,
    double U[], double R[], double V[], double T[],
    double E[], double P[], double Q[], double H[], double C_ad[])
{

    double time_current;
    int i;

    int inner_zones = zones / 10;
    int outer_zones = zones - inner_zones;

    double R_blast = .01 * pc;

    printf("Initializing blast at innermost %d cells (R = %le pc). \n", inner_zones, R_blast / pc);
    // INITIALIZE DATA FOR A NEW RUN
    for (i=0; i<zones; i++)
    {
        // set defaults
        U[i] = 0;

        // R[i] = i * (R_total / zones);

        double M = M_SN;
        // M = (1./V_0) * (4.*M_PI/3.) * pow(R_blast,3);
        
        if (i < inner_zones)
        {
            R[i] = i * R_blast / inner_zones;
            V[i] = (4.*M_PI/3)*pow(R_blast,3.) / M;
            E[i] = E_SN * M_SN / M;
            // V[i] = V_0 + ( ((4.*M_PI/3)*pow(R_blast,3.) / M) - V_0 ) * (inner_zones - i) / inner_zones;
            // E[i] = (c_V * T_0) + ( E_SN - (c_V * T_0)) * (inner_zones - i) / inner_zones;
            T[i] = E[i] / c_V;
            P[i] = (gamma - 1) * E[i] / V[i];
            Q[i] = 0;
            H[i] = 0;
            C_ad[i] = sqrt(gamma * P[i] * V[i]);
        }
        else
        {
            R[i] = R_blast
                + (i-inner_zones) * (R_total - R_blast) / (zones - inner_zones);
            V[i] = V_0;
            T[i] = T_0;
            E[i] = c_V * T[i];
            P[i] = (gamma - 1) * E[i] / V[i];
            Q[i] = 0;
            H[i] = 0;
            C_ad[i] = sqrt(gamma * P[i] * V[i]);
        }
    }

    time_current = 0;
    return time_current;
}

double rt1d_comparison(int zones,
    double U[], double R[], double V[], double T[],
    double E[], double P[], double Q[], double H[], double C_ad[])
{

    double time_current;
    int i;

    double r1 = 1e-5;
    double r2 = .5;

    double log_r1 = log(r1);
    double log_r2 = log(r2);

    double rSN = 1e-2;

    printf("Setting up initial conditions identical to RT1D sedov run \n");

    // INITIALIZE DATA FOR A NEW RUN
    for (i=0; i<zones; i++)
    {
        // set defaults
        U[i] = 0;

        R[i] = r1 * exp( (i * (log_r2 - log_r1) / zones) );
        
        if (R[i] < rSN)
        {
            V[i] = 1e-2;
            T[i] = 1e5;
        }
        else
        {
            V[i] = 1.0;
            T[i] = 1.0;
        }
        
        E[i] = c_V * T[i];
        P[i] = (gamma - 1) * E[i] / V[i];
        Q[i] = 0;
        H[i] = 0;
        C_ad[i] = sqrt(gamma * P[i] * V[i]);
    }


    time_current = 0;
    return time_current;
}


double restart(int *zones_ptr, int *k_ptr,
        double U[], double R[], double V[], double T[],
        double E[], double P[], double Q[], double H[], double C_ad[],
        char mod_filename[])
{
    double time_current; 


    FILE *mod_file_pointer; // for restarting runs
          mod_file_pointer    = fopen(mod_filename, "r"); 



    // READ DATA TO CONTINUE PREVIOUS RUN
    char buffer[1024];  // buffer for reading file
    printf("\nReading in from file: %s \n", mod_filename);

    fgets(buffer, sizeof(buffer), mod_file_pointer); // skip headers
    fscanf(mod_file_pointer,"%lf \n", &time_current);
    fgets(buffer, sizeof(buffer), mod_file_pointer); // skip headers
    fscanf(mod_file_pointer,"%d \n", k_ptr);
    fgets(buffer, sizeof(buffer), mod_file_pointer); // skip headers
    fscanf(mod_file_pointer,"%d \n", zones_ptr);
    fgets(buffer, sizeof(buffer), mod_file_pointer); // skip headers


    printf("starting time: %lf [yr]\n", time_current / yr);
    int zones = *zones_ptr;
    int i;
    for (i=0; i<zones; ++i)
    {
        fscanf(mod_file_pointer, 
            "%lf %lf %lf %lf %lf %lf %lf %lf %lf \n", 
            &(R[i]),
            &(U[i]),
            &V[i],
            &T[i],
            &E[i],
            &P[i],
            &Q[i],
            &H[i],
            &C_ad[i]); 

    }
    fclose(mod_file_pointer);


    return time_current;
}