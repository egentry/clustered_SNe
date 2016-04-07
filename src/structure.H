#ifndef _STRUCTURE_H
#define _STRUCTURE_H

// respectively: density, Pressure, radial velocity, metallicity
enum{RHO,PPP,VRR,ZZZ};

// respectively: mass, energy, momentum (radial)
enum{DDD,TAU,SRR}; 

#include <string>

#include <stdlib.h>
#include <string.h>
#include <time.h>
extern "C" {
#include <grackle.h> // cooling_units
}

#include <vector>

#define NUM_Q 4 // number of evolved variables
#define NUM_G 1 // number of guard cells per end

struct param_list{

    double  t_min, t_max;

    int     NumRepts, NumChecks;
    bool    Out_LogTime;

    int     Num_R;
    double  rmin, rmax;
    int     LogZoning;
    double  LogRadius;

    double  MaxShort, MaxLong;

    double  CFL;
    int     PLM; // Piece-wise reconstruction: 1 = piecewise linear, 0 = piecewise constant
    int     RK2; // Integration: 1 = 2nd order Runge-Kutta; 0 = 1st order Eulerian
    double  H_0, H_1; // artificial conduction constants

    int     Riemann_Solver;
    int     Mesh_Motion;
    double  Density_Floor, Pressure_Floor;

    bool    With_Cooling;
    std::string cooling_type;
    double  Cooling_Redshift;
    
    double  Adiabatic_Index;

    std::string ICs;
    
    std::string mass_loss;

};

struct domain{

    struct cell * theCells;
    int Nr,Ng;                 // number of cells (total and guard, respectively)

    time_t Wallt_init;
    std::string output_prefix;

    struct param_list theParList;  // for reading in from "in.par" file -- see above

    double t;
    int count_steps;
    double t_init, t_fin;

    // n*  = currently at the nth action of *
    // N_* = user desires N actions of * 
    // (some may be skipped if the timestep jumps over multiple * actions)
    int nrpt;
    int N_rpt;
    int nsnp;    // unused
    int N_snp;   // unused
    int nchk;
    int N_chk;
    int nchk_0; // checkpoint numbering starts *at* nchk_0 (used for restart)

    int final_step;        // flag that lets the code know when to terminate

    // For grackle cooling: 
    code_units cooling_units;    

    // For initial conditions:
    double metallicity;              
    double background_density;       
    double background_temperature;   

    std::vector<struct supernova>  SNe;

    double cluster_mass; // [solar masses]
    int seed;

    double R_shock;


};

struct cell{
    double prim[NUM_Q];   // primitive    variables  (density, pressure,       velocity, etc.)
    double cons[NUM_Q];   // conservative variables  (Mass,    Energy (total), momentum, etc.)
    double RKcons[NUM_Q]; // holds the value of the conservative variables before a Runge-kutta timestep is started
    double grad[NUM_Q];   // spatial gradient, for PLM reconstruction
    double riph;          // radius of outside edge of cell
    double dr;            // total dr of cell (difference between outside edge + inside edge of cell)
    double wiph;          // velocity of right edge of cell

    // for use in checking dE_adiabatic (in fix_negative_energies):
    double E_int_old; 
    double dV_old; 
    double dE_cool; // negative for cooling, positive for heating
};


struct supernova
{
    double mass;           // [g]
    double mass_ejecta;    // [g]
    double mass_ejecta_Z;  // [g]
    double mass_winds;     // [g]
    double lifetime;       // [s]
};


#endif
