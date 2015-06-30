
#include <string.h>
#include <assert.h>
#include <math.h>
#include <grackle.h>

#include "structure.h"

double get_dA( double );
double get_dV( double , double );
double get_moment_arm( double , double );

double mindt( double * , double , double , double );

double getmindt( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   double dt = 1e100;
   int i;
   for( i=1 ; i<Nr-1 ; ++i ){
      int im = i-1;
      struct cell * c = theCells+i;
      double dr = c->dr;
      double r = c->riph-.5*dr;
      double wm = theCells[im].wiph;
      double wp = c->wiph;
      double w = .5*(wm+wp);
      double dt_temp = mindt( c->prim , w , r , dr );
      if( dt > dt_temp ) dt = dt_temp;
   }
   dt *= theDomain->theParList.CFL; 
   MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );

   return( dt );
}

void prim2cons( double * , double * , double );
void cons2prim( double * , double * , double );
double get_vr( double * );

void set_wcell( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int mesh_motion = theDomain->theParList.Mesh_Motion;
   int Nr = theDomain->Nr;

   int i;
   for( i=0 ; i<Nr-1 ; ++i ){
      struct cell * cL = theCells+i;  
      double w = 0.0;
      if( mesh_motion ){
         struct cell * cR = theCells+i+1;
         double wL = get_vr( cL->prim );
         double wR = get_vr( cR->prim );
         w = .5*(wL + wR); 
         if( i==0 && theDomain->rank==0 ) w = wR*(cR->riph - .5*cR->dr)/(cL->riph);//0.0;//2./3.*wR;
      }
      cL->wiph = w;
   }
}

void adjust_RK_cons( struct domain * theDomain , double RK ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
      }

      double prim_tmp[NUM_Q];
      double rp = c->riph;
      double rm = rp-c->dr;
      double dV = get_dV( rp , rm );
      cons2prim( c->cons , prim_tmp , dV );

      if(prim_tmp[PPP] < theDomain->theParList.Pressure_Floor)
      {
         printf("------ ERROR in adjust_RK_cons()------- \n");
         printf("pressure should be above pressure floor! \n");
         printf("pressure       = %e in cell %i \n", prim_tmp[PPP], i);
         printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
         printf("rp = %e \n", rp);
         printf("rm = %e \n", rm);
         printf("dr = %e \n", c->dr);
         assert(0);
      }
      // int q;
      for( q=0 ; q<NUM_Q ; ++q)
      {
         if(!isfinite(prim_tmp[q]) && q!=AAA)
         {
            printf("------ ERROR in adjust_RK_cons()------- \n");
            printf("prim[%d] = %e in cell %d \n", q, prim_tmp[q], i);
            printf("rp = %e \n", rp);
            printf("rm = %e \n", rm);
            printf("dr = %e \n", c->dr);
            assert(0);
         }
      }

   }

}

void move_cells( struct domain * theDomain , double RK , double dt){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;

      double prim_tmp[NUM_Q];
      double rp = c->riph;
      double rm = rp-c->dr;
      double dV = get_dV( rp , rm );
      cons2prim( c->cons , prim_tmp , dV );
      if(prim_tmp[PPP] < theDomain->theParList.Pressure_Floor)
      {
         printf("------ ERROR in move_cells() before moving------- \n");
         printf("pressure should be above pressure floor! \n");
         printf("pressure       = %e in cell %i \n", prim_tmp[PPP], i);
         printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
         printf("rp = %e \n", c->riph);
         printf("rm = %e \n", c->riph - c->dr);
         printf("dr = %e \n", c->dr);
         assert(0);
      }

      c->riph += c->wiph*dt;

      // verify postconditions
      rp = c->riph;
      rm = rp-c->dr;
      dV = get_dV( rp , rm );
      cons2prim( c->cons , prim_tmp , dV );
      if(prim_tmp[PPP] < theDomain->theParList.Pressure_Floor)
      {
         printf("------ ERROR in move_cells() after moving------- \n");
         printf("pressure should be above pressure floor! \n");
         printf("pressure       = %e in cell %i \n", prim_tmp[PPP], i);
         printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
         printf("rp = %e \n", rp);
         printf("rm = %e \n", rm);
         printf("dr = %e \n", c->dr);
         assert(0);
      }
      int q;
      for( q=0 ; q<NUM_Q ; ++q)
      {
         if(!isfinite(prim_tmp[q]) && q!=AAA)
         {
            printf("------ ERROR in move_cells()------- \n");
            printf("non-finite prim[q] \n");
            printf("prim[%d] = %e in cell %d \n", q, prim_tmp[q], i);
            printf("rp = %e \n", rp);
            printf("rm = %e \n", rm);
            printf("dr = %e \n", c->dr);
            assert(0);
         }
      }

   }

}

void calc_dr( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;
   for( i=1 ; i<Nr ; ++i ){
      int im = i-1;
      double rm = theCells[im].riph;
      double rp = theCells[i ].riph;
      double dr = rp-rm;
      theCells[i].dr = dr;
   }
   if( theDomain->rank==0 ) theCells[0].dr = theCells[0].riph;

}

void fix_negative_energies( struct domain * theDomain )
{
   // Only TOTAL energy is conserved
   //    - this can lead to a problem if kinetic energy exceeds total energy
   //    - in order to compensate, internal energy (and pressure) will become negative
   // To fix this problem:
   //    - we specifically evolve internal energy using adiabatic assumptions
   //    - if the adiabatic internal energy differs from the numeric internal energy by too much
   //       then we use the adiabatic energy rather than the numeric energy

   double tolerance = .5; // relative tolerance allowed on energy error before switching to adiabatic

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double gamma = theDomain->theParList.Adiabatic_Index;

   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;

      // double E_old = c->P_old * c->dV_old / (gamma - 1);

      double rp = c->riph;
      double rm = rp-c->dr;
      double dV = get_dV( rp , rm );

      double P_new = c->P_old * pow(dV / c->dV_old, -1*gamma);
      double E_new =    P_new * dV / (gamma - 1);

      double Mass = c->cons[DDD];
      double vr   = c->cons[SRR] / Mass;
      double E_numeric = c->cons[TAU] - .5*Mass*vr*vr;

      if ( (((E_new - E_numeric) / E_new) > tolerance) || (E_numeric < 0) ) 
      {
         // overwrite the energy, so that it has a strictly positive internal energy
         printf("------ Applying energy fix for cell %d  ------- \n", i);
         printf("time = %e \n", theDomain->t);
         printf("P_old = %e \n", c->P_old);
         printf("dV_old = %e \n", c->dV_old);
         printf("dV     = %e \n", dV);
         printf("P_new = %e \n", P_new);
         printf("E_new = %e \n", E_new);
         printf("P_numeric = %e \n", E_numeric * (gamma-1) / dV);
         printf("E_numeric = %e \n", E_numeric);
         c->cons[TAU] = E_new + .5*Mass*vr*vr;

         c->P_old  = P_new;
         c->dV_old = dV;


      } 
      else
      {
         c->P_old  = E_numeric * (gamma - 1) / dV;
         c->dV_old = dV;
      }


   }

   return;
}


void calc_prim( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double dV = get_dV( rp , rm );
      cons2prim( c->cons , c->prim , dV );

      // verify postconditions
      if(c->prim[PPP] < theDomain->theParList.Pressure_Floor)
      {
         printf("------ ERROR in calc_prim()------- \n");
         printf("pressure should be above pressure floor! \n");
         printf("pressure       = %e in cell %i \n", c->prim[PPP], i);
         printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
         printf("rp = %e \n", rp);
         printf("rm = %e \n", rm);
         printf("dr = %e \n", c->dr);
         assert(0);
      }
      int q;
      for( q=0 ; q<NUM_Q ; ++q)
      {
         if(!isfinite(c->prim[q]) && q!=AAA)
         {
            printf("------ ERROR in calc_prim()------- \n");
            printf("prim[%d] = %e in cell %d \n", q, c->prim[q], i);
            printf("rp = %e \n", rp);
            printf("rm = %e \n", rm);
            printf("dr = %e \n", c->dr);
            assert(0);
         }
      }
   }
}

void plm( struct domain *);
void riemann( struct cell * , struct cell * , double , double , double );

void radial_flux( struct domain * theDomain , double dt ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   plm( theDomain );
   for( i=0 ; i<Nr-1 ; ++i ){
      struct cell * cL = theCells+i;
      struct cell * cR = theCells+i+1;
      double r = cL->riph;
      double dA = get_dA(r); 
      riemann( cL , cR , r , dA , dt );
   }

}

void source( double * , double * , double , double , double , double , double , code_units , int );
void source_alpha( double * , double * , double * , double , double );

void add_source( struct domain * theDomain , double dt ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double grad[NUM_Q];

   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double r = get_moment_arm(rp,rm);
      double dV = get_dV(rp,rm);
      source( c->prim , c->cons , rp , rm , dV , dt , 
         theDomain->metallicity ,  theDomain->cooling_units , theDomain->theParList.With_Cooling );
      int inside = i>0 && i<Nr-1;
      for( q=0 ; q<NUM_Q ; ++q ){
         if( inside ){
            struct cell * cp = theCells+i+1;
            struct cell * cm = theCells+i-1;
            double dR = .5*cp->dr + c->dr + .5*cm->dr;
            grad[q] = (cp->prim[q]-cm->prim[q])/dR;
         }else{
            grad[q] = 0.0;
         }
      }
      source_alpha( c->prim , c->cons , grad , r , dV*dt );


      // verify postconditions
      if(c->prim[PPP] < theDomain->theParList.Pressure_Floor)
      {
         printf("------ ERROR in add_source()------- \n");
         printf("pressure should be above pressure floor! \n");
         printf("pressure       = %e in cell %i \n", c->prim[PPP], i);
         printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
         printf("rp = %e \n", rp);
         printf("rm = %e \n", rm);
         printf("dr = %e \n", c->dr);
         assert(0);
      }
      for( q=0 ; q<NUM_Q ; ++q)
      {
         if(!isfinite(c->prim[q]) && q!=AAA)
         {
            printf("------ ERROR in add_source()------- \n");
            printf("prim[%d] = %e in cell %d \n", q, c->prim[q], i);
            printf("rp = %e \n", rp);
            printf("rm = %e \n", rm);
            printf("dr = %e \n", c->dr);
            assert(0);
         }
      }
   }   

}


void longandshort( struct domain * theDomain , double * L , double * S , int * iL , int * iS , int * rL , int * rS ){ 

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double rmax = theCells[Nr-1].riph;
   double rmin = theCells[0].riph;
   MPI_Allreduce( MPI_IN_PLACE , &rmax , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &rmin , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );
   int Nr0 = theDomain->theParList.Num_R;
   double dr0 = rmax/(double)Nr0;
   double dx0 = log(rmax/rmin)/Nr0;
   int logscale = theDomain->theParList.LogZoning;

   double Long  = 0.0; 
   double Short = 0.0; 
   int iLong  = -1;
   int iShort = -1;

   int rank = theDomain->rank;
   int size = theDomain->size;
   int Ng   = theDomain->Ng;

   int imin = 1;
   int imax = Nr-1;
   if( rank!=0 )      imin = Ng;
   if( rank!=size-1 ) imax = Nr-Ng;

   int i;
   for( i=imin ; i<imax ; ++i ){
      struct cell * c = theCells+i;
      double dy = c->dr;
      double dx = dr0;
      if( logscale ) dx = c->riph*dx0;
      double l = dy/dx;
      double s = dx/dy;
      if( Long  < l ){ Long  = l; iLong  = i; } 
      if( Short < s ){ Short = s; iShort = i; } 
   }

   struct { double value ; int index ; } maxminbuf;
   maxminbuf.value = Short;
   maxminbuf.index = rank;
   MPI_Allreduce( MPI_IN_PLACE , &maxminbuf , 1 , MPI_DOUBLE_INT , MPI_MAXLOC , MPI_COMM_WORLD );
   *rS = maxminbuf.index;

   maxminbuf.value = Long;
   maxminbuf.index = rank;
   MPI_Allreduce( MPI_IN_PLACE , &maxminbuf , 1 , MPI_DOUBLE_INT , MPI_MAXLOC , MPI_COMM_WORLD );
   *rL = maxminbuf.index;

   *iS = iShort;
   *iL = iLong;
   *S = Short;
   *L = Long;

}


void AMR( struct domain * theDomain ){

   double L,S;
   int iL=0;
   int iS=0;
   int rL=0;
   int rS=0;
   longandshort( theDomain , &L , &S , &iL , &iS , &rL , &rS );
   int rank = theDomain->rank;

   //if( rank == rL ) printf("Rank %d; Long  = %e #%d\n",rank,L,iL);
   //if( rank == rS ) printf("Rank %d; Short = %e #%d\n",rank,S,iS);

   double MaxShort = theDomain->theParList.MaxShort;
   double MaxLong  = theDomain->theParList.MaxLong;

   struct cell * theCells = theDomain->theCells;
   int Ng = theDomain->Ng;
   int Nr = theDomain->Nr;

   if( S>MaxShort && rank == rS ){
      // printf("KILL!  iS = %d\n",iS);

      int iSp = iS+1;
      int iSm = iS-1;
      //Possibly shift iS backwards by 1 
      double drL = theCells[iSm].dr;
      double drR = theCells[iSp].dr;
      int imin = Ng;
      if( rank==0 ) imin = 0;
      if( drL<drR && iSm>imin ){
         --iS;
         --iSm;
         --iSp;
      }
      struct cell * c  = theCells+iS;
      struct cell * cp = theCells+iSp;

      //Remove Zone at iS+1
      c->dr   += cp->dr;
      c->riph  = cp->riph;
      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q]   += cp->cons[q];
         c->RKcons[q] += cp->RKcons[q];
      }
      double rp = c->riph;
      double rm = rp - c->dr;
      double dV = get_dV( rp , rm );
      cons2prim( c->cons , c->prim , dV );
      c->P_old  = c->prim[PPP];
      c->dV_old = dV;
      //Shift Memory
      int blocksize = Nr-iSp-1;
      memmove( theCells+iSp , theCells+iSp+1 , blocksize*sizeof(struct cell) );
      theDomain->Nr -= 1;
      Nr = theDomain->Nr;
      theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
      theCells = theDomain->theCells;
      if( iS < iL ) iL--;

   }

   if( L>MaxLong && rank==rL ){

      // printf("FORGE! iL = %d\n",iL);
      theDomain->Nr += 1;
      Nr = theDomain->Nr;
      theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
      theCells = theDomain->theCells;
      int blocksize = Nr-iL-1;
      memmove( theCells+iL+1 , theCells+iL , blocksize*sizeof(struct cell) );

      struct cell * c  = theCells+iL;
      struct cell * cp = theCells+iL+1;

      double rp = c->riph;
      double rm = rp - c->dr;
      double r0 = pow( .5*(rp*rp*rp+rm*rm*rm) , 1./3. );

      c->riph  = r0;
      c->dr    = r0-rm;
      cp->dr   = rp-r0;

      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q]    *= .5;
         c->RKcons[q]  *= .5;
         cp->cons[q]   *= .5;
         cp->RKcons[q] *= .5;
      }

      double dV = get_dV( r0 , rm );
      cons2prim( c->cons , c->prim , dV );
      c->P_old = c->prim[PPP];
      c->dV_old = dV;
      dV = get_dV( rp , r0 );
      cons2prim( cp->cons , cp->prim , dV );
      cp->P_old = cp->prim[PPP];
      cp->dV_old = dV;

   }

}



