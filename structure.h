#ifndef _STRUCTURE_H
#define _STRUCTURE_H

enum{RHO,PPP,VRR,XXX,AAA}; // respectively: density, density * Pressure, radial velocity, XXX (passive scalar), AAA (eddy strength)
enum{DDD,TAU,SRR}; // respectively: mass, energy, momentum (radial)?

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <grackle.h>

#define NUM_Q 5
#define NUM_G 2

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

   double CFL, PLM;
   double Density_Floor, Pressure_Floor;
   int With_Cooling;
   double Adiabatic_Index;

};

struct domain{

   struct cell * theCells;
   int Nr,Ng;

   time_t Wallt_init;
   int rank,size;

   struct param_list theParList;

   double t;
   int count_steps;
   double t_init, t_fin;
   int nrpt;
   int N_rpt;
   int nsnp;
   int N_snp;
   int nchk;
   int N_chk;

   int final_step;

   code_units cooling_units;
   double metallicity;

};

struct cell{
   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double grad[NUM_Q];
   double riph;
   double dr;
   double wiph;
   // for use in checking dE_adiabatic
   double P_old; 
   double dV_old; 
};


#endif
