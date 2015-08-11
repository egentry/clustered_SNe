
#include <assert.h>
#include <time.h>

#include <grackle.h>

#include "../structure.h"

double get_moment_arm( double , double );
double get_dV( double , double );

void output( struct domain * theDomain , const char * filestart , double t ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   char filename[256] = "";
   strcat(filename, theDomain->output_prefix);
   strcat(filename, filestart);
   strcat(filename, ".dat");

   FILE * pFile = fopen( filename , "w" );
   fprintf(pFile,"# time = %le [s] \n", t);
   // fprintf(pFile,"# r           dr           dV           Density      Pressure     Velocity     X            Alpha\n");
   fprintf(pFile,"# r                  dr                 dV                 Density            Pressure           Velocity           X                  Alpha\n");

   int i_min = 0;
   int i_max = Nr;

   int i,q;
   for( i=i_min ; i<i_max ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double dr = c->dr; 
      double rm = rp-dr;
      double dV = get_dV( rp , rm );
      fprintf(pFile,"%18.10e %18.10e %18.10e ",rp,dr,dV);
      for( q=0 ; q<NUM_Q ; ++q ){
         fprintf(pFile,"%18.10e ",c->prim[q]);
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

void overview( struct domain * theDomain )
{
   // prints an overview of key parameters into a datafile
   
   char overview_filename[256] = "";
   strcat(overview_filename, theDomain->output_prefix);
   strcat(overview_filename, "overview.dat");
   FILE * oFile = fopen(overview_filename,"w");

   fprintf(oFile, "Metallicity:            %e \n", theDomain->metallicity);
   fprintf(oFile, "Background Density:     %e \n", theDomain->background_density);
   fprintf(oFile, "Background Temperature: %e \n", theDomain->background_temperature);
   fprintf(oFile, "With cooling:           %d \n", theDomain->theParList.With_Cooling);
   time_t current_time = time(NULL);
   fprintf(oFile, "Created at: %s \n", ctime(&current_time));

   fclose(oFile);
}
