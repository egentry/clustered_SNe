
#ifndef _READPAR_H_
#define _READPAR_H_

int readvar( std::string filename , const char * varname , int vartype , 
		     void * ptr );

int read_par_file( struct domain * theDomain , int argc , char * argv [] );



#endif