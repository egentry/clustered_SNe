
#include <iostream>
#include <string.h>
#include <string>

#include "structure.H"
#include "readpar.H"

enum{VAR_INT,VAR_DOUB,VAR_STR};

int readvar( std::string filename , const char * varname , int vartype , 
             void * ptr ){

   FILE * inFile = fopen( filename.c_str() , "r" );
   if (inFile==NULL)
   {
      std::cerr << "Input file \"" << filename << "\" not found." << std::endl;
      return 1;
   }

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
   if( found==0 ){
       std::cerr << "Variable not found by readvar: " << varname << std::endl;
      return 1;
   } 

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

   return 0;
}

int read_par_file( struct domain * theDomain , int argc , char * argv [] ){

    struct param_list * theList = &( theDomain->theParList );

    std::string par_filename("in.par");
    if ( argc > 1 )
    {
        par_filename = std::string(argv[1]);
    }

    char tmp_str[512];

    int err=0;  

    err += readvar( par_filename , "Num_R"           , VAR_INT  ,
                                  &(theList->Num_R)           );
    err += readvar( par_filename , "Num_Reports"     , VAR_INT  ,
                                  &(theList->NumRepts)        );
    err += readvar( par_filename , "Num_Snapshots"   , VAR_INT  ,
                                  &(theList->NumSnaps)        );
    err += readvar( par_filename , "Num_Checkpoints" , VAR_INT  ,
                                  &(theList->NumChecks)       );
    err += readvar( par_filename , "T_Start"         , VAR_DOUB ,
                                  &(theList->t_min)           );
    err += readvar( par_filename , "T_End"           , VAR_DOUB ,
                                  &(theList->t_max)           );
    err += readvar( par_filename , "R_Min"           , VAR_DOUB ,
                                  &(theList->rmin)            );
    err += readvar( par_filename , "R_Max"           , VAR_DOUB ,
                                  &(theList->rmax)            );
    err += readvar( par_filename , "Use_Logtime"     , VAR_INT  ,
                                  &(theList->Out_LogTime)     );
    err += readvar( par_filename , "Log_Zoning"      , VAR_INT  ,
                                  &(theList->LogZoning)       );
    err += readvar( par_filename , "Log_Radius"      , VAR_DOUB ,
                                  &(theList->LogRadius)       );
    err += readvar( par_filename , "CFL"             , VAR_DOUB ,
                                  &(theList->CFL)             );
    err += readvar( par_filename , "PLM"             , VAR_INT  ,
                                  &(theList->PLM)             );
    err += readvar( par_filename , "RK2"             , VAR_INT  ,
                                  &(theList->RK2)             );
    err += readvar( par_filename , "Adiabatic_Index" , VAR_DOUB ,
                                  &(theList->Adiabatic_Index) );
    err += readvar( par_filename , "Density_Floor"   , VAR_DOUB ,
                                  &(theList->Density_Floor)   );
    err += readvar( par_filename , "Pressure_Floor"  , VAR_DOUB ,
                                  &(theList->Pressure_Floor)  );
    err += readvar( par_filename , "With_Cooling"    , VAR_INT  ,
                                  &(theList->With_Cooling)    );
    err += readvar( par_filename , "Mesh_Motion"     , VAR_INT  ,
                                  &(theList->Mesh_Motion)     );
    err += readvar( par_filename , "Riemann_Solver"  , VAR_INT  ,
                                  &(theList->Riemann_Solver)  );
    err += readvar( par_filename , "Use_RT"          , VAR_INT  ,
                                  &(theList->rt_flag)         );
    err += readvar( par_filename , "Use_Logtime"     , VAR_INT  ,
                                  &(theList->Out_LogTime)     );
    err += readvar( par_filename , "Max_Aspect_Short", VAR_DOUB ,
                                  &(theList->MaxShort)        );
    err += readvar( par_filename , "Max_Aspect_Long" , VAR_DOUB ,
                                  &(theList->MaxLong)         );
    err += readvar( par_filename , "ICs" , VAR_STR ,
                                  tmp_str                     );

    theList->ICs = std::string(tmp_str);

    if( err > 0 )
    {
        std::cerr << "Reading parameter file failed" << std::endl;
        return 1;
    }

   return 0;

}


