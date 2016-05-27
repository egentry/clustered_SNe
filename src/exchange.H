#ifndef _exchange_H_
#define _exchange_H_

#include <mpi.h>

void generate_mpi_cell( MPI_Datatype * cell_mpi );
void copy_cell_to_lite( struct cell * c , struct cell_lite * cl );
void copy_lite_to_cell( struct cell_lite * cl , struct cell * c ); 
void generate_sendbuffer( struct domain * theDomain , struct cell_lite * pl , 
    struct cell_lite * pr , int mode );
void exchange_data( struct domain * theDomain );


#endif