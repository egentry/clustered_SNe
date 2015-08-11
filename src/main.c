
#include "structure.h"

int    read_par_file(            struct domain * );
int    parse_command_line_args ( struct domain * , int , char *[] );
void   setupGrid(                struct domain * );
int    setupDomain(              struct domain * );
void   setupCells(               struct domain * );
void   overview(                 struct domain * );
void   start_clock(              struct domain * );
void   set_wcell(                struct domain * );
double getmindt(                 struct domain * );
void   check_dt(                 struct domain * , double * );
void   possiblyOutput(           struct domain * , int );
void   timestep(                 struct domain * , double );
void   generate_log(             struct domain * );
void   freeDomain(               struct domain * );

int main( int argc , char * argv[] ){
 
   int error;

   struct domain theDomain = {0};

   error = read_par_file( &theDomain );
   if( error==1 ) 
   {
      return(0);
   }
   
   error =  parse_command_line_args ( &theDomain , argc , argv );
   if( error==1 ) 
   {
      return(0);
   }

   setupGrid( &theDomain );   
   error = setupDomain( &theDomain );
   if( error==1 ) 
   {
      return(0);
   }

   setupCells( &theDomain );

   overview( &theDomain );

   start_clock( &theDomain ); 

   while( !(theDomain.final_step) ){

      set_wcell( &theDomain );
      double dt = getmindt( &theDomain );
      check_dt( &theDomain , &dt );
      possiblyOutput( &theDomain , 0 );
      timestep( &theDomain , dt );

   }

   possiblyOutput( &theDomain , 1 );

   generate_log( &theDomain );
   freeDomain( &theDomain );

   return(0);

}