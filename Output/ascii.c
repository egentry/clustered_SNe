#include <assert.h>
#include "../paul.h"

double get_moment_arm( double , double );
double get_dV( double , double );

void output( struct domain * theDomain , char * filestart , double t){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   int rank = theDomain->rank;
   int size = theDomain->size;

   char filename[256];
   sprintf(filename,"%s.dat",filestart);

   if( rank==0 ){
      FILE * pFile = fopen( filename , "w" );
      fprintf(pFile,"# time = %le [s] \n", t);
      fprintf(pFile,"# r           dr           dV           Density      Pressure     Velocity     X            Alpha\n");
      fclose(pFile);
   }
   MPI_Barrier( MPI_COMM_WORLD );

   int i_min = 0;
   int i_max = Nr;

   if( rank != 0      ) i_min = Ng;
   if( rank != size-1 ) i_max = Nr-Ng;

   int rk;
   for( rk=0 ; rk<size ; ++rk ){
      if( rank==rk ){
         FILE * pFile = fopen( filename , "a" );
         int i,q;
         for( i=i_min ; i<i_max ; ++i ){
            struct cell * c = theCells+i;
            double rp = c->riph;
            double dr = c->dr; 
            double rm = rp-dr;
            double r  = get_moment_arm( rp , rm );
            double dV = get_dV( rp , rm );
            fprintf(pFile,"%e %e %e ",r,dr,dV);
            for( q=0 ; q<NUM_Q ; ++q ){
               fprintf(pFile,"%e ",c->prim[q]);
            }
            if(c->prim[PPP] < theDomain->theParList.Pressure_Floor)
            {
               printf("-------ERROR--------- \n");
               printf("In output() \n");
               printf("c->prim[PPP] = %e at i=%d \n", c->prim[PPP], i);
               printf("Pressure floor should be = %e \n", theDomain->theParList.Pressure_Floor);
               assert(0);
            }
            fprintf(pFile,"\n");
         }
         fclose( pFile );
      }
      MPI_Barrier( MPI_COMM_WORLD );
   }

}
