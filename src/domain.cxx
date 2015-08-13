
#include <iostream>

#include <uuid/uuid.h>
#include <grackle.h>

#include <assert.h>

#include "domain.H"
#include "structure.H"
#include "constants.H"
#include "blast.H"

#include <algorithm>

code_units setup_cooling( struct domain * );

int  setICparams( struct domain * );
void setHydroParams( struct domain * );
void setRiemannParams( struct domain * );

int setupDomain( struct domain * theDomain ){

   int error;

   if (theDomain->theParList.With_Cooling == 1)
   {
      theDomain->cooling_units = setup_cooling(theDomain);
   }


   theDomain->t       = theDomain->theParList.t_min;
   theDomain->t_init  = theDomain->theParList.t_min;
   theDomain->t_fin   = theDomain->theParList.t_max;

   theDomain->N_rpt = theDomain->theParList.NumRepts;
   theDomain->N_snp = theDomain->theParList.NumSnaps;
   theDomain->N_chk = theDomain->theParList.NumChecks;

   theDomain->count_steps = 0;
   theDomain->final_step  = 0;

   theDomain->nrpt=-1;
   theDomain->nsnp=-1;
   theDomain->nchk=-1;

   // // if you don't want a uuid prefix: 
   // strcpy(theDomain->output_prefix, "test_");
   
   uuid_t  id_binary;
   char id_ascii[80];
   uuid_generate(id_binary);
   uuid_unparse(id_binary, id_ascii);
   printf("uuid: %s \n", id_ascii);
   strcat(id_ascii, "_");
   strcpy(theDomain->output_prefix, id_ascii);

   error = setICparams( theDomain );
   if ( error==1 ) return(error);
   setHydroParams( theDomain );
   setRiemannParams( theDomain );

   // everything below should probably go in a separate function?

   assert( std::is_sorted( theDomain->SNe_times.rbegin(),
                           theDomain->SNe_times.rend()   ) );

   // SNe_times should have been set within the initial conditions
   // probably within parse_command_line_args?
   if ( theDomain->SNe_times.size() > 0 )
   {

      double t_first_SN = theDomain->SNe_times.back();
      double t_last_SN  = theDomain->SNe_times.front();

      theDomain->t      += t_first_SN;
      theDomain->t_init += t_first_SN;
      theDomain->t_fin  += t_last_SN;

   }
   else
   {
      std::cerr << "Error: No SNe in this run. Exiting." << std::endl;
      // no supernovae. For now, just kill the process
      // but maybe figure out a better way to respond?
      return 1; 
   }

   return(0);

}

double get_moment_arm( double , double );
double get_dV( double , double );
void initial( double * , double ); 
void prim2cons( double * , double * , double );
void cons2prim( double * , double * , double );
void boundary( struct domain * );
void setupCells( struct domain * theDomain ){

   int i;
   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = &(theCells[i]);
      double rp = c->riph;
      double rm = rp - c->dr;
      c->wiph = 0.0; 
      double r = get_moment_arm( rp , rm );
      double dV = get_dV( rp , rm );
      initial( c->prim , r ); 
      prim2cons( c->prim , c->cons , dV );
      cons2prim( c->cons , c->prim , dV );
      c->P_old  = c->prim[PPP];
      c->dV_old = dV;
   }

   boundary( theDomain );

}

void freeDomain( struct domain * theDomain ){
   free( theDomain->theCells );
}

void check_dt( struct domain * theDomain , double * dt ){

   double t = theDomain->t;
   double tmax = theDomain->t_fin;
   int final=0;
   if( t + *dt > tmax ){
      *dt = tmax-t;
      final=1;
   }

   FILE * abort = NULL;
   abort = fopen("abort","r");
   if( abort ){ final = 1; fclose(abort); }

   if( final ) theDomain->final_step = 1;

}

// void snapshot( struct domain * , char * );
void output( struct domain * , const char *, double );

void possiblyOutput( struct domain * theDomain , int override ){

   double t = theDomain->t;
   double t_min = theDomain->t_init;
   double t_fin = theDomain->t_fin;
   double Nrpt = theDomain->N_rpt;
   // double Nsnp = theDomain->N_snp;
   double Nchk = theDomain->N_chk;
   int LogOut = theDomain->theParList.Out_LogTime;
   int n0;

   n0 = (int)( t*Nrpt/t_fin );
   if( LogOut ) n0 = (int)( Nrpt*log(t/t_min)/log(t_fin/t_min) );
   if( theDomain->nrpt < n0 || override ){
      theDomain->nrpt = n0;
      printf("t = %.3e\n",t);
   }

   n0 = (int)( t*Nchk/t_fin );
   if( LogOut ) n0 = (int)( Nchk*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nchk < n0 && Nchk>0) || override ){
      theDomain->nchk = n0;
      char filename[256];
      if( !override ){
         printf("Creating Checkpoint #%04d...\n",n0);
         sprintf(filename,"checkpoint_%04d",n0);
         output( theDomain , filename, t );
      }else{
         printf("Creating Final Checkpoint...\n");
         output( theDomain , "output", t );
      }
   }
/*
   n0 = (int)( t*Nsnp/t_fin );
   if( LogOut ) n0 = (int)( Nsnp*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nsnp < n0 && Nsnp>0) || override ){
      theDomain->nsnp = n0;
      char filename[256];
      if(!override) sprintf( filename , "snapshot_%04d" , n0 );
      else sprintf( filename , "snapshot" );
      // snapshot( theDomain , filename );
   }
*/
}


