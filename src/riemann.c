
#include <math.h>
#include <assert.h>

#include "structure.h"

enum{_HLL_,_HLLC_};

int riemann_solver = 0;
int rt_flag = 0;

static double PRE_FLOOR = 0.0;


void setRiemannParams( struct domain * theDomain ){
   riemann_solver = theDomain->theParList.Riemann_Solver;
   rt_flag = theDomain->theParList.rt_flag;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
}

void prim2cons( double * , double * , double );
void cons2prim( double * , double * , double );
void flux( double * , double * );
void getUstar( double * , double * , double , double );
void vel( double * , double * , double * , double * , double * , double );
double get_eta( double * , double * , double );
double get_dV( double , double );


void riemann( struct cell * cL , struct cell * cR, double r , double dA , double dt ){

   double dAdt = dA * dt;
   double primL[NUM_Q];
   double primR[NUM_Q];

   double drL = .5*cL->dr;
   double drR = .5*cR->dr;

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + cL->grad[q]*drL;
      primR[q] = cR->prim[q] - cR->grad[q]*drR;
   }

   // ======== Verify pre-conditions ========= //
   if (primL[PPP] < PRE_FLOOR)
   {
      printf("------ERROR in riemann()------- \n");
      printf("primL[%d] = %e \n", PPP, primL[PPP]);
      printf("Floor should be: %e \n", PRE_FLOOR);
      printf("cL->prim[%d] = %e \n", PPP, cL->prim[PPP]);
      printf("Occured at r = %e \n", r);
      printf("drL = %e \n", drL);
      printf("drR = %e \n", drR);
      printf("cL->grad[PPP] = %e \n", cL->grad[PPP]);
   }

   if (primR[PPP] < PRE_FLOOR)
   {
      printf("------ERROR in riemann() ------- \n");
      printf("primR[%d] = %e \n", PPP, primR[PPP]);
      printf("Floor should be: %e \n", PRE_FLOOR);
      printf("cR->prim[%d] = %e \n", PPP, cR->prim[PPP]);
      printf("Occured at r = %e \n", r);
      printf("drL = %e \n", drL);
      printf("drR = %e \n", drR);
      printf("cR->grad[PPP] = %e \n", cR->grad[PPP]);

      assert(0);
   }

   // continue with Riemann solver

   double Sl,Sr,Ss;

   vel( primL , primR , &Sl , &Sr , &Ss , r );

   double Fl[NUM_Q];
   double Fr[NUM_Q];
   double Ul[NUM_Q];
   double Ur[NUM_Q];

   double Flux[NUM_Q];

   double w = cL->wiph;

   if( w < Sl ){
      flux( primL , Fl );
      prim2cons( primL , Ul , 1.0 );  // using dV = 1 gives the cons per unit volume

      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fl[q] - w*Ul[q];
      if(!isfinite(Flux[q]) && q!=AAA)
            {
               printf("------ERROR in riemann()------- \n");
               printf("Bad flux in part 0 of riemann() \n");
            }
      }
   }else if( w > Sr ){
      flux( primR , Fr );
      prim2cons( primR , Ur , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fr[q] - w*Ur[q];
         if(!isfinite(Flux[q]) && q!=AAA)
         {
            printf("------ERROR in riemann()------- \n");
            printf("Bad flux in part 1 of riemann() \n");
         }
      }
   }else{
      if( riemann_solver == _HLL_ ){
         double Fstar;
         double Ustar;
         double aL =  Sr;
         double aR = -Sl;
 
         prim2cons( primL , Ul , 1.0 );
         prim2cons( primR , Ur , 1.0 );
         flux( primL , Fl );
         flux( primR , Fr );

         for( q=0 ; q<NUM_Q ; ++q ){
            Fstar = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
            Ustar = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );

            Flux[q] = Fstar - w*Ustar;
            if(!isfinite(Flux[q]) && q!=AAA)
            {
               printf("------ERROR in riemann()------- \n");
               printf("Bad flux in part 2 of riemann() \n");
            }
         }
      }else{
         double Ustar[NUM_Q];
         double Uk[NUM_Q];
         double Fk[NUM_Q];
         if( w < Ss ){
            prim2cons( primL , Uk , 1.0 );
            getUstar( primL , Ustar , Sl , Ss ); 
            flux( primL , Fk ); 

            for( q=0 ; q<NUM_Q ; ++q ){
               Flux[q] = Fk[q] + Sl*( Ustar[q] - Uk[q] ) - w*Ustar[q];
               if(!isfinite(Flux[q]) && q!=AAA)
               {
                  printf("------ERROR in riemann()------- \n");
                  printf("Bad flux in part 3 of riemann() \n");
               }
            }    
         }else{
            prim2cons( primR , Uk , 1.0 );
            getUstar( primR , Ustar , Sr , Ss ); 
            flux( primR , Fk ); 

            for( q=0 ; q<NUM_Q ; ++q ){
               Flux[q] = Fk[q] + Sr*( Ustar[q] - Uk[q] ) - w*Ustar[q];
               if(!isfinite(Flux[q]) && q!=AAA)
               {
                  printf("------ERROR in riemann()------- \n");
                  printf("Bad flux in part 4, Flux[%d]=%e \n", q, Flux[q]);
                  printf("Fk[%d] = %e \n", q, Fk[q]);
                  printf("Sr    = %e \n", Sr);
                  printf("Ustar[%d] = %e \n", q, Ustar[q]);
                  printf("w = %e \n", w);
               }
            } 
         } 
      }
   }

   for ( q=0 ; q<NUM_Q ; ++q)
   {
      if(!isfinite(Flux[q]) && q!=AAA)
      {
         printf("------ERROR in riemann()------- \n");
         printf("Non-finite Flux[%d]= %e at r = %e \n", q, Flux[q], r);
         printf("cL->dr = %e \n", cL->dr);
         printf("cR->dr = %e \n", cR->dr);
         assert(0);
      }
   }

   if( rt_flag ){
      double prim[NUM_Q];
      double consL[NUM_Q];
      double consR[NUM_Q];
      prim2cons( cL->prim , consL , 1.0 );
      prim2cons( cR->prim , consR , 1.0 );
      double gprim[NUM_Q];
      double gcons[NUM_Q];
      for( q=0 ; q<NUM_Q ; ++q ){
         prim[q] = .5*(primL[q]+primR[q]);
         gprim[q] = (cR->prim[q] - cL->prim[q])/(drL+drR);
         gcons[q] = (consR[q] - consL[q])/(drL+drR);
      }
      double eta = get_eta( prim , gprim , r );
      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] += -eta*gcons[q];
      }
   }



   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dAdt;
      cR->cons[q] += Flux[q]*dAdt;
   }


   // ======== Verify post-conditions ========= //
   double primL_tmp[NUM_Q];
   double primR_tmp[NUM_Q];

   double rp_L = cL->riph;
   double rm_L = rp_L - cL->dr;
   double dV_L = get_dV( rp_L , rm_L );
   cons2prim( cL->cons , primL_tmp , dV_L );  // using dV = 1 gives the cons per unit volume -- counteracts dV = 1 above
   double rp_R = cR->riph;
   double rm_R = rp_R - cR->dr;
   double dV_R = get_dV( rp_R , rm_R );
   cons2prim( cR->cons , primR_tmp , dV_R);
   for( q=0 ; q<NUM_Q ; ++q)
   {
      if(!isfinite(primL_tmp[q]) && q!=AAA)
      {
         printf("------ERROR in riemann()------- \n");
         printf("primL[%d] not finite by end of riemann() \n", q);
         printf("primL[%d] = %e \n", q, primL_tmp[q]);
         assert(0);
      }
      if(!isfinite(primR_tmp[q]) && q!=AAA)
      {
         printf("------ERROR in riemann()------- \n");
         printf("primR[%d] not finite by end of riemann() \n", q);
         printf("primR[%d] = %e \n", q, primR_tmp[q]);
         assert(0);
      }
   }
   if(primL_tmp[PPP] < PRE_FLOOR)
   {
      printf("------ERROR in riemann()------- \n");
      printf("while preparing to exit \n");
      printf("left prim[%d] = %e \n", PPP, primL_tmp[PPP]);
      printf("expected P > PRE_FLOOR \n");
      assert(0);
   }
   if(primR_tmp[PPP] < PRE_FLOOR)
   {
      printf("------ERROR in riemann()------- \n");
      printf("while preparing to exit \n");
      printf("right prim[%d] = %e \n", PPP, primR_tmp[PPP]);
      printf("expected P > PRE_FLOOR \n");
      assert(0);
   }

}


