#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <constants.h>

#include <grid.h>


void swap(double **a, double **b)
{
    // must be dynamically allocated array

    double *tmp = *a;

    *a = *b;
    *b = tmp;

    return;
}

void expand_grid(int *zones_ptr, 
        double **U_old_ptr, double **U_new_ptr,
        double **R_old_ptr, double **R_new_ptr,
        double **V_old_ptr, double **V_new_ptr,
        double **T_old_ptr, double **T_new_ptr,
        double **E_old_ptr, double **E_new_ptr, double **dE_ptr,
        double **P_old_ptr, double **P_new_ptr, 
        double **Q_old_ptr, double **Q_new_ptr,
        double **H_old_ptr, double **H_new_ptr,
        double **C_ad_ptr,
        double **Mass_ptr,  double **M_int_ptr)
{
    

    double *U_old_tmp, *U_new_tmp;
    double *R_old_tmp, *R_new_tmp;
    double *V_old_tmp, *V_new_tmp;
    double *T_old_tmp, *T_new_tmp;
    double *E_old_tmp, *E_new_tmp, *dE_tmp;
    double *P_old_tmp, *P_new_tmp;
    double *Q_old_tmp, *Q_new_tmp;
    double *H_old_tmp, *H_new_tmp;
    double *C_ad_tmp;
    double *Mass_tmp, *M_int_tmp;

    int zones     = *(zones_ptr);

    int new_zones = (zones * 2) - 1; // MUST STAY ODD


    U_old_tmp = malloc(new_zones * sizeof(double));
    U_new_tmp = malloc(new_zones * sizeof(double));
    R_old_tmp = malloc(new_zones * sizeof(double));
    R_new_tmp = malloc(new_zones * sizeof(double));
    V_old_tmp = malloc(new_zones * sizeof(double));
    V_new_tmp = malloc(new_zones * sizeof(double));
    T_old_tmp = malloc(new_zones * sizeof(double));
    T_new_tmp = malloc(new_zones * sizeof(double));
    E_old_tmp = malloc(new_zones * sizeof(double));
    E_new_tmp = malloc(new_zones * sizeof(double));
    dE_tmp    = malloc(new_zones * sizeof(double));
    P_old_tmp = malloc(new_zones * sizeof(double));
    P_new_tmp = malloc(new_zones * sizeof(double));
    Q_old_tmp = malloc(new_zones * sizeof(double));
    Q_new_tmp = malloc(new_zones * sizeof(double));
    H_old_tmp = malloc(new_zones * sizeof(double));
    H_new_tmp = malloc(new_zones * sizeof(double));
    C_ad_tmp  = malloc(new_zones * sizeof(double));
    Mass_tmp  = malloc(new_zones * sizeof(double));
    M_int_tmp = malloc(new_zones * sizeof(double));


    double dR = ((*R_old_ptr)[zones-2] - (*R_old_ptr)[zones-4]) / 2;

    int i;
    for (i=0; i<new_zones; ++i)
    {

        if (i < zones-2)
        {
            // If within the old range, just copy the old value
            U_old_tmp[i]    = (*U_old_ptr)[i];
            U_new_tmp[i]    = (*U_new_ptr)[i];
            R_old_tmp[i]    = (*R_old_ptr)[i];
            R_new_tmp[i]    = (*R_new_ptr)[i];
            V_old_tmp[i]    = (*V_old_ptr)[i];
            V_new_tmp[i]    = (*V_new_ptr)[i];
            T_old_tmp[i]    = (*T_old_ptr)[i];
            T_new_tmp[i]    = (*T_new_ptr)[i];
            E_old_tmp[i]    = (*E_old_ptr)[i];
            E_new_tmp[i]    = (*E_new_ptr)[i];
            dE_tmp   [i]    = (*dE_ptr)   [i];
            P_old_tmp[i]    = (*P_old_ptr)[i];
            P_new_tmp[i]    = (*P_new_ptr)[i];
            Q_old_tmp[i]    = (*Q_old_ptr)[i];
            Q_new_tmp[i]    = (*Q_new_ptr)[i];
            H_old_tmp[i]    = (*H_old_ptr)[i];
            H_new_tmp[i]    = (*H_new_ptr)[i];            
            C_ad_tmp [i]    = (*C_ad_ptr) [i];
        }
        else
        {
            // If we need to make a new cell, just extrapolate outer boundary
            U_old_tmp[i]    = 0;
            U_new_tmp[i]    = 0;
            R_old_tmp[i]    = (*R_old_ptr)[zones-4] + (i-(zones-4))*dR;
            R_new_tmp[i]    = (*R_new_ptr)[zones-4] + (i-(zones-4))*dR;
            V_old_tmp[i]    = (*V_old_ptr)[zones-3];
            V_new_tmp[i]    = (*V_new_ptr)[zones-3];
            T_old_tmp[i]    = (*T_old_ptr)[zones-3];
            T_new_tmp[i]    = (*T_new_ptr)[zones-3];
            E_old_tmp[i]    = (*E_old_ptr)[zones-3];
            E_new_tmp[i]    = (*E_new_ptr)[zones-3];
            dE_tmp   [i]    = 0;
            P_old_tmp[i]    = (*P_old_ptr)[zones-3];
            P_new_tmp[i]    = (*P_new_ptr)[zones-3];
            Q_old_tmp[i]    = 0;
            Q_new_tmp[i]    = 0;
            H_old_tmp[i]    = 0;
            H_new_tmp[i]    = 0;            
            C_ad_tmp [i]    = (*C_ad_ptr) [zones-3];
        }
    }   

    (*zones_ptr) = new_zones;
    zones = *zones_ptr;

    // Overwrite Mass, M_int without checking
    create_Mass_array(Mass_tmp, R_old_tmp, V_old_tmp, zones);
    create_M_int_array(M_int_tmp, Mass_tmp, zones); 

    swap(&U_old_tmp, U_old_ptr);
    swap(&U_new_tmp, U_new_ptr);
    swap(&R_old_tmp, R_old_ptr);
    swap(&R_new_tmp, R_new_ptr);
    swap(&V_old_tmp, V_old_ptr);
    swap(&V_new_tmp, V_new_ptr);
    swap(&T_old_tmp, T_old_ptr);
    swap(&T_new_tmp, T_new_ptr);
    swap(&E_old_tmp, E_old_ptr);
    swap(&E_new_tmp, E_new_ptr);
    swap(&dE_tmp,    dE_ptr);
    swap(&P_old_tmp, P_old_ptr);
    swap(&P_new_tmp, P_new_ptr);
    swap(&Q_old_tmp, Q_old_ptr);
    swap(&Q_new_tmp, Q_new_ptr);
    swap(&H_old_tmp, H_old_ptr);
    swap(&H_new_tmp, H_new_ptr);
    swap(&C_ad_tmp,  C_ad_ptr);
    swap(&Mass_tmp,  Mass_ptr);
    swap(&M_int_tmp, M_int_ptr);


    free(U_old_tmp);
    free(U_new_tmp);
    free(R_old_tmp);
    free(R_new_tmp);
    free(V_old_tmp);
    free(V_new_tmp);
    free(T_old_tmp);
    free(T_new_tmp);
    free(E_old_tmp);
    free(E_new_tmp);
    free(dE_tmp);
    free(P_old_tmp);
    free(P_new_tmp);
    free(Q_old_tmp);
    free(Q_new_tmp);
    free(H_old_tmp);
    free(H_new_tmp);
    free(C_ad_tmp);
    free(Mass_tmp);
    free(M_int_tmp);

    return;


}


void create_Mass_array(double Mass[], double R[], double V[], int zones)
{
    int i;
    for(i=0; i<zones-2; i+=2)
    {
        if (i==0)
        {
            Mass[i] = (4./3) * M_PI * pow(R[i+1],3) / V[i];
        }
        else
        {
            Mass[i] = (4./3) * M_PI * (pow(R[i+1],3) - pow(R[i-1],3))
                / V[i];
        }
    }


    for(i=1; i<zones-1; i+=2)
    {
        Mass[i] = interp(Mass[i-1], Mass[i+1]);
    }
 
    return;
};

void create_M_int_array(double M_int[], double Mass[], int zones)
{

    int i;
    //  Calculate Mass internal to each edge
    M_int[1] = Mass[0];
    for(i=3; i<zones-1; i+=2)
    {
        M_int[i] = M_int[i-2] + Mass[i-1];
    }    

    return;
}

void enforce_boundary_conditions_Neumann(int zones, 
    double P[], 
    double Q[] )
{
    //  IMPLICIT
    //      Dirichlet:
    //          R - inside and outside
    //          U - inside and outside
    //          E - inside and outside
    //          T - 

    //  EXPLICIT
    //      Neumann:
    //          - Pressure      - inside and outside
    //          - Viscosity     - inside and outside

    //  NONE
    //      Energy (internal) - inside and outside
    //      Temperature - inside and outside

    // U should probably be in Dirichlet function,
    //      but it needs to be done for for the ICs as well

    // U[0] = 0; // should be implicitly enforced by neumann on P,Q 
    // U[1] = 0; 
    P[0] = P[2];
    // P[1] = P[2];
    Q[0] = Q[2];
    // Q[1] = Q[2];


    // U[zones-1] = 0; // should be implicitly enforced by neumann on P,Q 
    // U[zones-2] = 0; 
    P[zones-1] = P[zones-3];
    // P[zones-2] = P[zones-3];
    Q[zones-1] = Q[zones-3];
    // Q[zones-2] = Q[zones-3];


    return;
}

void enforce_boundary_conditions_Dirichlet(int zones, 
    double U_old[], double U_new[], 
    double R_old[], double R_new[], 
    double V_old[], double V_new[], 
    double T_old[], double T_new[], 
    double E_old[], double E_new[],
    double H_old[], double H_new[] )
{
    //  IMPLICIT
    //      Dirichlet:
    //          R - inside and outside
    //          U - inside and outside
    //          E - inside and outside
    //          T - 

    //  EXPLICIT
    //      Neumann:
    //          - Pressure      - inside and outside
    //          - Viscosity     - inside and outside

    //  NONE
    //      Energy (internal) - inside and outside
    //      Temperature - inside and outside

    // REMEMBER, if you need to change U,
    //      you need to change it in the *_Neumann function as well

    // U_new[0] = 0; // should be implicitly enforced by neumann on P,Q 
    // U_new[1] = 0; 


    // V_new[0] = V_old[0];
    // T_new[0] = T_old[0];
    // E_new[0] = E_old[0];

    // R_new[1] = R_old[1];


    // U_new[zones-1] = 0; // should be implicitly enforced by neumann on P,Q 
    // V_new[zones-1] = V_old[zones-1];
    // T_new[zones-1] = T_old[zones-1];
    // E_new[zones-1] = E_old[zones-1];

    // U_new[zones-2] = 0; 
    // R_new[zones-2] = R_old[zones-2];


    H_old[1] = 0;
    H_new[1] = 0;
    H_old[zones-2] = 0;
    H_new[zones-2] = 0;


    return;
}

void enforce_boundary_conditions(int zones, 
    double U_old[], double U_new[], 
    double R_old[], double R_new[], 
    double V_old[], double V_new[], 
    double T_old[], double T_new[], 
    double E_old[], double E_new[], 
    double P_old[], double P_new[], 
    double Q_old[], double Q_new[],
    double H_old[], double H_new[] )
{
    //  IMPLICIT
    //      Dirichlet:
    //          R - inside and outside
    //          U - inside and outside
    //          E - inside and outside
    //          T - 

    //  EXPLICIT
    //      Neumann:
    //          - Pressure      - inside and outside
    //          - Viscosity     - inside and outside

    //  NONE
    //      Energy (internal) - inside and outside
    //      Temperature - inside and outside

    enforce_boundary_conditions_Neumann(zones, 
        P_new, 
        Q_new );

    enforce_boundary_conditions_Dirichlet(zones, 
        U_old, U_new, 
        R_old, R_new, 
        V_old, V_new, 
        T_old, T_new, 
        E_old, E_new,
        H_old, H_new );
}

void enforce_boundary_conditions_initial(int zones, 
    double U_old[], double U_new[], 
    double R_old[], double R_new[], 
    double V_old[], double V_new[], 
    double T_old[], double T_new[], 
    double E_old[], double E_new[], 
    double P_old[], double P_new[], 
    double Q_old[], double Q_new[],
    double H_old[], double H_new[] )
{
    //  IMPLICIT
    //      Dirichlet:
    //          R - inside and outside
    //          U - inside and outside
    //          E - inside and outside
    //          T - 

    //  EXPLICIT
    //      Neumann:
    //          - Pressure      - inside and outside
    //          - Viscosity     - inside and outside

    //  NONE
    //      Energy (internal) - inside and outside
    //      Temperature - inside and outside

    U_old[1]        = 0;
    U_old[zones-2]  = 0;
    U_new[1]        = 0;
    U_new[zones-2]  = 0;

    Q_old[0]        = 0; // shouldn't be necessary
    Q_old[zones-1]  = 0; // shouldn't be necessary
    Q_new[0]        = 0; // shouldn't be necessary
    Q_new[zones-1]  = 0; // shouldn't be necessary

    H_old[1]        = 0;
    H_old[zones-2]  = 0;
    H_new[1]        = 0;
    H_new[zones-2]  = 0;

    R_new[0]        = R_old[0];
    R_new[1]        = R_old[1];
    R_new[zones-2]  = R_old[zones-2];
    R_new[zones-1]  = R_old[zones-1];

    V_new[0]        = V_old[0];
    V_new[1]        = V_old[1];
    V_new[zones-2]  = V_old[zones-2];
    V_new[zones-1]  = V_old[zones-1];



}


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

double harmonic_mean(double a, double b)
{
    double mean;
 
    mean = 2 * (a*b)/(a+b);

    return mean;
}