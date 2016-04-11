
#include <assert.h>
#include <iostream>
#include <string.h>
#include <string>

#include "structure.H"
#include "readpar.H"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>

int read_var_raw( const std::string filename , const char * var_name ,  char * var_raw )
{

    FILE * inFile = fopen( filename.c_str() , "r" );
    if (inFile==NULL)
    {
        std::cerr << "Input file \"" << filename << "\" not found." << std::endl;
        return 1;
    }

    char line[512];
    bool found = false;
    std::vector<std::string> strs;

    while( (fgets(line,512,inFile) != NULL) && found==false )
    {

        boost::split( strs , line , boost::is_any_of("\t\n "), 
                      boost::token_compress_on);
        if( strcmp(strs[0].c_str() , var_name) == 0 )
        {
            strcpy( var_raw , strs[1].c_str() );
            found = true;
        }
    }

    fclose( inFile );
 
    if( found == false )
    {
        std::cerr << "Variable not found by read_var_raw: " << var_name << std::endl;

        strcpy( var_raw , "" );
        return 1;  
    } 

    return 0;
}

int read_var( const std::string filename , const char * var_name , 
              int * var , const int var_default )
{

    char var_raw[512];
    int error = read_var_raw( filename , var_name , var_raw );

    if ( error == 0 )
    {
        *var = std::stoi(std::string(var_raw));
    }
    else
    {
        std::cerr << "For variable: "   << var_name << " "
                  << "using default = " << var_default << std::endl;

        *var = var_default; 
    }

    return 0;

}

int read_var( const std::string filename , const char * var_name , 
              bool * var , const bool var_default )
{

    char var_raw[512];
    int error = read_var_raw( filename , var_name , var_raw );
    std::string var_raw_string = var_raw;
    boost::to_lower(var_raw_string);

    if ( error == 0 )
    {
        if (var_raw_string.compare(std::string("true")) == 0)
        {
            *var = true;
        }
        else if (var_raw_string.compare("false") == 0 )
        {
            *var = false;
        }
        else
        {
            *var = boost::lexical_cast<bool>(var_raw);
        }
    }
    else
    {
        *var = var_default; 
    }

    return 0;

}

int read_var( const std::string filename , const char * var_name , 
              double * var , const double var_default )
{

    char var_raw[512];
    int error = read_var_raw( filename , var_name , var_raw );

    if ( error == 0 )
    {
        *var = std::stod(std::string(var_raw));
    }
    else
    {
        std::cerr << "For variable: "   << var_name << " "
                  << "using default = " << var_default << std::endl;

        *var = var_default; 
    }

    return 0;

}

int read_var( const std::string filename , const char * var_name , 
              std::string * var , const std::string var_default )
{

    char var_raw[512];
    int error = read_var_raw( filename , var_name , var_raw );

    if ( error == 0 )
    {
        *var = std::string(var_raw);
    }
    else
    {
        std::cerr << "For variable: "   << var_name << " "
                  << "using default = " << var_default << std::endl;

        *var = var_default; 
    }

    return 0;

}

int read_par_file( struct domain * theDomain , int argc , char * argv [] )
{

    struct param_list * theList = &( theDomain->theParList );

    std::string par_filename("in.par");
    if ( argc > 1 )
    {
        par_filename = std::string(argv[1]);
    }

    int error = 0;  

    error += read_var( par_filename , 
                    "T_Start"          , &(theList->t_min) , 0.0 );
    error += read_var( par_filename , 
                    "T_End"            , &(theList->t_max) , 1.0 );
    error += read_var( par_filename , 
                    "Num_Reports"      , &(theList->NumRepts) , 1000 );
    error += read_var( par_filename , 
                    "Num_Checkpoints"  , &(theList->NumChecks) , 100 );
    error += read_var( par_filename , 
                    "Use_Logtime"      , &(theList->Out_LogTime) , true );
    error += read_var( par_filename , 
                    "Num_R"            , &(theList->Num_R) , 1024 );
    error += read_var( par_filename , 
                    "R_Min"            , &(theList->rmin) , 3.08e16 );
    error += read_var( par_filename , 
                    "R_Max"            , &(theList->rmax) , 3.08e20 );
    error += read_var( par_filename , 
                    "Log_Zoning"       , &(theList->LogZoning) , 0 );
    error += read_var( par_filename , 
                    "Log_Radius"       , &(theList->LogRadius) , 3.08e18 );
    error += read_var( par_filename , 
                    "Max_Aspect_Short" , &(theList->MaxShort) , 10.0 );
    error += read_var( par_filename , 
                    "Max_Aspect_Long"  , &(theList->MaxLong) , 10.0 );
    error += read_var( par_filename , 
                    "CFL"              , &(theList->CFL) , 0.2 );
    error += read_var( par_filename , 
                    "PLM"              , &(theList->PLM) , 1 );
    error += read_var( par_filename , 
                    "RK2"              , &(theList->RK2) , 1 );
    error += read_var( par_filename , 
                    "H_0"              , &(theList->H_0) , 0.0 );
    error += read_var( par_filename , 
                    "H_1"              , &(theList->H_1) , 1.0 );
    error += read_var( par_filename , 
                    "Riemann_Solver"   , &(theList->Riemann_Solver) , 1 );
    error += read_var( par_filename , 
                    "Mesh_Motion"      , &(theList->Mesh_Motion) , 1 );
    error += read_var( par_filename , 
                    "Density_Floor"    , &(theList->Density_Floor) , 1e-60 );
    error += read_var( par_filename , 
                    "Pressure_Floor"   , &(theList->Pressure_Floor) , 1e-40 );
    error += read_var( par_filename , 
                    "With_Cooling"     , &(theList->With_Cooling) , 1 );
    error += read_var( par_filename , 
                    "Cooling_Type"     , &(theList->cooling_type) , std::string("equilibrium") );
    error += read_var( par_filename , 
                    "Cooling_Redshift" , &(theList->Cooling_Redshift) , 0.0 );
    error += read_var( par_filename , 
                    "Adiabatic_Index"  , &(theList->Adiabatic_Index) , 1.666666667 );
    error += read_var( par_filename , 
                    "ICs"              , &(theList->ICs) , std::string("cluster_SNe") );
    error += read_var( par_filename , 
                    "mass_loss"        , &(theList->mass_loss) , std::string("uniform") );


    if( error > 0 )
    {
        std::cerr << "Reading parameter file failed" << std::endl;
        return 1;
    }

   return 0;

}


