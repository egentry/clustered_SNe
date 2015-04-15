
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <constants.h> // defines physical constants

#include <ICs.h>


double sedov(int zones, double *U, double *R, double *V, double *T, 
            double *E, double *P, double *Q, double *C_ad)
{
    double time_current; 

    int sedov_zones = zones / 5;
    int outer_zones = zones - sedov_zones;
    char command[80];
    int ret;
    ret = sprintf(command, "python dimensionalize_sedov_initial.py %d %d %le", sedov_zones, outer_zones, R_total / 20 );
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

    printf("starting time: %lf [s]\n", time_current);
    int i;
    for (i=0; i<zones; ++i)
    {
        fscanf(in_file_pointer, 
            "%lf %lf %lf %lf %lf %lf %lf %lf \n", 
            &(U[i]),
            &(R[i]),
            &V[i],
            &T[i],
            &E[i],
            &P[i],
            &Q[i],
            &C_ad[i]); 

    }
    fclose(in_file_pointer);

    return time_current;
}


void new_blast(int zones, int zones_inner,
    double U[], double R[], double V[], double T[],
    double E[], double P[], double Q[], double C_ad[])
{
    int i;

    // INITIALIZE DATA FOR A NEW RUN
    for (i=0; i<zones; i++)
    {
        // set defaults
        U[i] = 0;  // this might be overwritten for inner zones
        Q[i] = 0;

        if (i==0)
        {
            //  UNRESOLVED INNER CORE (zone)
            //      effectively acts as a guard cell

            V[i] = V_core;
            T[i] = T_core;
            E[i] = E_core; 
            P[i] = 0;           // get from neumann boundary conditions
            Q[i] = Q_core;
            C_ad[ i] = C_ad_core;
        }
        else if (i==1)
        {
            //  UNRESOLVED INNER CORE (edge)
            //      effectively acts as a guard cell
            R[i] = R_core;

            // these are a little poorly defined for edges
            //      but it sets overall boundary conditions
            V[i] = V_core;
            T[i] = T_core;
            E[i] = E_core; 
            P[i] = 0;           // get from neumann boundary conditions
            Q[i] = Q_core;
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
                    T[i] = T_inner;

                }
                else
                {
                    //  EDGES
                    R[i] = i / ((double) zones_inner) * R_inner;

                }
                // variables which don't particularly care 
                //  if they're in the zone or edge
                V[i] = V_inner;
                T[i] = T_inner;
                E[i] = E_SN;
                P[i] = P_inner;
                C_ad[ i] = sqrt(gamma * P[i] * V[i]);
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
                    // background EDGES
                    R[i] = R_inner
                        + ( (i-zones_inner) / ((double) zones-zones_inner)
                            * (R_total-R_inner) );
                }

                // background variables which could be zones or edges
                V[i] = V_0;
                T[i] = T_0;
                E[i] = c_V * T_0;
                P[i] = comp * T_0 / V_0;
                C_ad[ i] = sqrt(gamma * P[i] * V[i]);
            }                
        }
    }

}

double restart(int *zones_ptr, int *k_ptr,
        double U[], double R[], double V[], double T[],
        double E[], double P[], double Q[], double C_ad[],
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


    printf("starting time: %lf [s]\n", time_current);
    int zones = *zones_ptr;
    int i;
    for (i=0; i<zones; ++i)
    {
        fscanf(mod_file_pointer, 
            "%lf %lf %lf %lf %lf %lf %lf %lf \n", 
            &(R[i]),
            &(U[i]),
            &V[i],
            &T[i],
            &E[i],
            &P[i],
            &Q[i],
            &C_ad[i]); 

    }
    fclose(mod_file_pointer);


    return time_current;
}