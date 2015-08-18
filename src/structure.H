#ifndef _STRUCTURE_H
#define _STRUCTURE_H

enum{RHO,PPP,VRR,XXX,AAA}; // respectively: density, Pressure, radial velocity, XXX (passive scalar), AAA (eddy strength)
enum{DDD,TAU,SRR}; // respectively: mass, energy, momentum (radial)

#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <grackle.h>

#include <vector>

#define NUM_Q 5 // number of evolved variables
#define NUM_G 2 // number of guard cells (counting both ends)

struct param_list{

   int Num_R;
   double t_min, t_max;
   double rmin,rmax;
   int NumRepts, NumSnaps, NumChecks;
   int Out_LogTime;

   int LogZoning;
   double LogRadius;
   int Mesh_Motion, Riemann_Solver;
   double MaxShort, MaxLong;
   int Absorb_BC, Initial_Regrid, rt_flag;

   double CFL;
   int    PLM; // Piece-wise reconstruction: 1 = piecewise linear, 0 = piecewise constant
   int    RK2; // Integration: 1 = 2nd order Runge-Kutta; 0 = 1st order Eulerian
   double Density_Floor, Pressure_Floor;
   int With_Cooling;
   double Adiabatic_Index;

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
   double metallicity;              

   // For Sedov initial conditions:
   double background_density;       
   double background_temperature;   

   std::vector<double> SNe_times;

   double cluster_mass; // [solar masses]
   int seed;


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
   double P_old; 
   double dV_old; 
};


#endif
