
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structure.H"


void calc_dr(struct domain * );
void setupGrid( struct domain * theDomain ){

   int Ng = NUM_G;
   theDomain->Ng = Ng;
   int Num_R = theDomain->theParList.Num_R;
   int LogZoning = theDomain->theParList.LogZoning;

   double Rmin = theDomain->theParList.rmin;
   double Rmax = theDomain->theParList.rmax;

   int Nr = Num_R;

   theDomain->Nr = Nr;
   theDomain->theCells = (struct cell *) malloc( Nr*sizeof(struct cell));

   int i;

   double dx = 1./(double)Num_R;
   double x0 = 0;
   double R0 = theDomain->theParList.LogRadius;
   for( i=0 ; i<Nr ; ++i ){
      double xp = x0 + ((double)i+1.)*dx;
      double rp;
      if( LogZoning == 0 ){
         rp = Rmin + xp*(Rmax-Rmin);
      }else if( LogZoning == 1 ){
         rp = Rmin*pow(Rmax/Rmin,xp);
      }else{
         rp = R0*pow(Rmax/R0,xp) + Rmin-R0 + (R0-Rmin)*xp;
      }
      theDomain->theCells[i].riph = rp;
   }
   calc_dr( theDomain );

}


