
#ifndef _READPAR_H_
#define _READPAR_H_

#include <string>

int read_var_raw( const std::string filename , const char * var_name ,  char * var_raw );

int read_var( const std::string filename , const char * var_name , 
              double * var , const double var_default = 0 );
int read_var( const std::string filename , const char * var_name , 
              int * var , const int var_default = 5 );
int read_var( const std::string filename , const char * var_name , 
              std::string * var , const std::string var_default = "" );
int read_var( const std::string filename , const char * var_name , 
              bool * var , const bool var_default );


int read_par_file( struct domain * theDomain , int argc , char * argv [] );



#endif