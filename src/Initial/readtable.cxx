
#include "../structure.H"

static int NL = 0;
static double * rr  = NULL;
static double * rho = NULL;
static double * Pp  = NULL;
static double * vr  = NULL;
static double * X   = NULL;
static double * A   = NULL;

int countlines(char * filename){
   FILE *pFile = fopen(filename, "r");
   int lines=0;
   char c;
   while ((c = fgetc(pFile)) != EOF){
      if (c == '\n') ++lines;
   }
   fclose(pFile);
   return(lines);
}

int getTable( void ){
   int nL = countlines("Initial/initial.dat") - 2;
   rr  = (double *) malloc( nL*sizeof(double) );
   double dr;
   double dV;
   rho = (double *) malloc( nL*sizeof(double) );
   Pp  = (double *) malloc( nL*sizeof(double) );
   vr  = (double *) malloc( nL*sizeof(double) );
   X   = (double *) malloc( nL*sizeof(double) );
   A   = (double *) malloc( nL*sizeof(double) );
   FILE * pFile = fopen("Initial/initial.dat","r");
   char tmp[1024];
   fgets(tmp, sizeof(tmp), pFile);
   // printf("%s\n", tmp);
   fgets(tmp, sizeof(tmp), pFile);
   // printf("%s\n", tmp);

   int l;
   for( l=0 ; l<nL ; ++l ){
      fscanf(pFile,"%le %le %le %le %le %le %le %le\n",&(rr[l]),&dr,&dV,&(rho[l]),&(Pp[l]),&(vr[l]),&(X[l]),&(A[l]));
      if(l==1)
      {
         printf("r = %le \n", rr[l]);
         printf("dr = %le \n", dr);
         printf("dV = %le \n", dV);
         printf("rho = %le \n", rho[l]);
         printf("P = %le \n", Pp[l]);
         printf("vr = %le \n", vr[l]);
         printf("X = %le \n", X[l]);         
         printf("A = %le \n", A[l]);         
      }

   }
   fclose(pFile);
   return(nL);
}

int setICparams( struct domain * theDomain ){
   NL = getTable();
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

int parse_command_line_args ( struct domain * theDomain , int argc , char * argv [] )
{
      return(0);
}