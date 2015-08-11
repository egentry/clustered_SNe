
#include <string.h>

#include "structure.h"

enum{VAR_INT,VAR_DOUB,VAR_STR};

int readvar( char * filename , const char * varname , int vartype , void * ptr ){

   FILE * inFile = fopen( filename , "r" );
   char s[512];
   char nm[512];
   char s1[512];
   int found = 0;
   
   while( (fgets(s,512,inFile) != NULL) && found==0 ){
      sscanf(s,"%s ",nm);
      if( strcmp(nm,varname)==0 ){
         strcpy(s1,s);
         found=1;
      }
   }
   
   fclose( inFile );
   if( found==0 ) return(1);

   char * s2 = s1+strlen(nm)+strspn(s1+strlen(nm),"\t :=>_");

   double temp;
   char stringval[256];

   sscanf(s2,"%lf",&temp);
   sscanf(s2,"%256s",stringval);

   if( vartype == VAR_INT ){
      *((int *)   ptr) = (int)temp;
   }else if( vartype == VAR_DOUB ){
      *((double *)ptr) = (double)temp;
   }else{
      strcpy( (char *) ptr , stringval );
   }

   return(0);
}

int read_par_file( struct domain * theDomain ){

   struct param_list * theList = &( theDomain->theParList );

   char pfile[] = "in.par";

   int err=0;  

   err += readvar( pfile , "Num_R"           , VAR_INT  , &(theList->Num_R)           );
   err += readvar( pfile , "Num_Reports"     , VAR_INT  , &(theList->NumRepts)        );
   err += readvar( pfile , "Num_Snapshots"   , VAR_INT  , &(theList->NumSnaps)        );
   err += readvar( pfile , "Num_Checkpoints" , VAR_INT  , &(theList->NumChecks)       );
   err += readvar( pfile , "T_Start"         , VAR_DOUB , &(theList->t_min)           );
   err += readvar( pfile , "T_End"           , VAR_DOUB , &(theList->t_max)           );
   err += readvar( pfile , "R_Min"           , VAR_DOUB , &(theList->rmin)            );
   err += readvar( pfile , "R_Max"           , VAR_DOUB , &(theList->rmax)            );
   err += readvar( pfile , "Use_Logtime"     , VAR_INT  , &(theList->Out_LogTime)     );
   err += readvar( pfile , "Log_Zoning"      , VAR_INT  , &(theList->LogZoning)       );
   err += readvar( pfile , "Log_Radius"      , VAR_DOUB , &(theList->LogRadius)       );
   err += readvar( pfile , "CFL"             , VAR_DOUB , &(theList->CFL)             );
   err += readvar( pfile , "PLM"             , VAR_INT  , &(theList->PLM)             );
   err += readvar( pfile , "RK2"             , VAR_INT  , &(theList->RK2)             );
   err += readvar( pfile , "Adiabatic_Index" , VAR_DOUB , &(theList->Adiabatic_Index) );
   err += readvar( pfile , "Density_Floor"   , VAR_DOUB , &(theList->Density_Floor)   );
   err += readvar( pfile , "Pressure_Floor"  , VAR_DOUB , &(theList->Pressure_Floor)  );
   err += readvar( pfile , "With_Cooling"    , VAR_INT  , &(theList->With_Cooling)    );
   err += readvar( pfile , "Mesh_Motion"     , VAR_INT  , &(theList->Mesh_Motion)     );
   err += readvar( pfile , "Riemann_Solver"  , VAR_INT  , &(theList->Riemann_Solver)  );
   err += readvar( pfile , "Use_RT"          , VAR_INT  , &(theList->rt_flag)         );
   err += readvar( pfile , "Use_Logtime"     , VAR_INT  , &(theList->Out_LogTime)     );
   err += readvar( pfile , "Max_Aspect_Short", VAR_DOUB , &(theList->MaxShort)        );
   err += readvar( pfile , "Max_Aspect_Long" , VAR_DOUB , &(theList->MaxLong)         );

   if( err > 0 ){
      printf("Read Failed\n");
      return(1);
   }

   return(0);

}


