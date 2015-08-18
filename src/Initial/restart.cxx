
#include <iostream>
#include <string>
#include <assert.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "../structure.H"

namespace fs = boost::filesystem;

static int NL = 0;
static double * rr  = NULL;
static double * rho = NULL;
static double * Pp  = NULL;
static double * vr  = NULL;
static double * X   = NULL;
static double * A   = NULL;

static std::string restart_filename("init");
static std::string restart_id;

static double time_restart = 0.0;
static unsigned int checkpoints_finished = 0;

int countlines(std::string filename){
   FILE *pFile = fopen(filename.c_str(), "r");
   if ( pFile == NULL )
   {
      std::cerr << "Error: Restart file (\"" << filename << "\" doesn't exist."
         << std::endl;
      return 0;
   }
   int lines=0;
   char c;
   while ((c = fgetc(pFile)) != EOF){
      if (c == '\n') ++lines;
   }
   fclose(pFile);
   return(lines);
}



int getTable( std::string filename ){
   int nL = countlines(filename) - 2;
   if ( nL < 0 ) return nL;
   rr  = (double *) malloc( nL*sizeof(double) );
   double dr;
   double dV;
   rho = (double *) malloc( nL*sizeof(double) );
   Pp  = (double *) malloc( nL*sizeof(double) );
   vr  = (double *) malloc( nL*sizeof(double) );
   X   = (double *) malloc( nL*sizeof(double) );
   A   = (double *) malloc( nL*sizeof(double) );
   FILE * pFile = fopen(filename.c_str(),"r");
   char tmp[1024];
   fscanf(pFile,"%s %s %s %le %s \n",
      tmp, tmp, tmp, &time_restart, tmp);
   std::cout << "restart time: " << time_restart << std::endl;
   
   fgets(tmp, sizeof(tmp), pFile); // header line

   int l;
   for( l=0 ; l<nL ; ++l ){
      fscanf(pFile,"%le %le %le %le %le %le %le %le\n",
         &(rr[l]),&dr,&dV,&(rho[l]),&(Pp[l]),&(vr[l]),&(X[l]),&(A[l]));
      // if(l==1)
      // {
      //    printf("r = %le \n", rr[l]);
      //    printf("dr = %le \n", dr);
      //    printf("dV = %le \n", dV);
      //    printf("rho = %le \n", rho[l]);
      //    printf("P = %le \n", Pp[l]);
      //    printf("vr = %le \n", vr[l]);
      //    printf("X = %le \n", X[l]);         
      //    printf("A = %le \n", A[l]);         
      // }

   }
   fclose(pFile);
   return(nL);
}

int setICparams( struct domain * theDomain ){
   NL = getTable( restart_filename ); 
   if ( NL < 0 )
   {
      std::cerr << "Error reading in restart file" << std::endl;
      return 1;
   }

   theDomain->t      += time_restart;
   theDomain->t_init += time_restart;
   theDomain->t_fin  += time_restart;
   theDomain->nchk_0 = checkpoints_finished + 1;

   if (restart_filename.find("checkpoint_") != std::string::npos )
   {
      const std::string restart_basename = fs::basename(restart_filename);
      const unsigned int id_start = restart_basename.find("checkpoint_");

      theDomain->output_prefix = restart_basename.substr(0,id_start);
   }




   return(0);
}

void initial( double * prim , double r ){

   // this routine (and the entire initial initial condition setup) is broken
   // I need a better way to be able to read in an old file into the domain.
   // Having initial conditions depend only on prim and r seems pretty broken

   int l=0;
   while( rr[l] < r && l < NL-2 ) ++l;
   if( l==0 ) ++l;

   double rp = rr[l];
   double rm = rr[l-1];
   double drm = fabs(r-rm);
   double drp = fabs(rp-r);

   double rh    = (rho[l-1]*drp + rho[l]*drm)/(drp+drm);
   double V     = (vr[l-1]*drp  + vr[l]*drm )/(drp+drm);
   double X = 1.0;

   if( l==NL-2 ){
      V = 0.0;
      X = 0.0;
      rh = 1./r/r;
   }
   double P = rh*1e-8;

   prim[RHO] = rh;
   prim[PPP] = P;
   prim[VRR] = V;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}

void freeTable( void ){
   free(rr);
   free(rho);
   free(Pp);
   free(vr);
   free(X);
   free(A);
}

std::string get_restart_filename( const std::string partial_ID )
{
   // there's probably a better way to do this.
   // boost::regex perhaps?

   int highest_checkpoint_seen = -1;
   
   fs::path current_dir("."); // make this an input parameter?

   fs::directory_iterator end_itr;

   // can I just auto this?

   for ( fs::directory_iterator dir_itr{current_dir} ;
         dir_itr != end_itr ; 
         ++dir_itr )
   {
      std::string tmp_filename = dir_itr->path().string();

      if (tmp_filename.find(partial_ID) == std::string::npos)
      {
         continue;
      }
      std::size_t pos_left = tmp_filename.find("checkpoint_");
      if (tmp_filename.find("checkpoint_") == std::string::npos )
      {
         continue;
      }
      pos_left += std::string("checkpoint_").size();

      std::size_t pos_right = tmp_filename.find(".dat");

      int checkpoint = std::stoi(tmp_filename.substr(pos_left, 
                                                     pos_right-pos_left));

      if ( checkpoint > highest_checkpoint_seen )
      {
         restart_filename.swap(tmp_filename);
         highest_checkpoint_seen = checkpoint;
      }
   }


   // check that file exists
   fs::path restart_path(restart_filename);
   assert(fs::exists(restart_path));

   std::cout << "restart filename: " << restart_filename << std::endl;

   if ( highest_checkpoint_seen >= 0 )
   {
      checkpoints_finished = highest_checkpoint_seen;
   }

   return restart_filename;

}

int parse_command_line_args ( struct domain * theDomain , int argc , 
                                                          char * argv [] )
{
   std::string partial_ID("A4");
   if ( argc > 2 )
   {
      partial_ID = std::string(argv[2]);
   }  

   restart_filename = get_restart_filename(partial_ID);

   return(0);
}