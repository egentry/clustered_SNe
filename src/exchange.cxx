
#include <mpi.h>

#include "structure.H"
#include "exchange.H"

struct cell_lite{
   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double grad[NUM_Q];   // spatial gradient, for PLM reconstruction
   double riph;
   double dr;
   double wiph;

    // for use in checking dE_adiabatic (in fix_negative_energies):
    double E_int_old; 
    double dV_old; 
    double dE_cool; // negative for cooling, positive for heating
};

void generate_mpi_cell( MPI_Datatype * cell_mpi ){

   struct cell_lite test;
   int count = 10;
   int blocksize[]      = {NUM_Q,NUM_Q,NUM_Q,NUM_Q,1,1,1,1,1,1};
   MPI_Datatype types[] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
   MPI_Aint offsets[10];

   offsets[0] = (char *)&(test.prim)   - (char *)(&test);
   offsets[1] = (char *)&(test.cons)   - (char *)(&test);
   offsets[2] = (char *)&(test.RKcons) - (char *)(&test);
   offsets[3] = (char *)&(test.grad)   - (char *)(&test);
   offsets[4] = (char *)&(test.riph)   - (char *)(&test);
   offsets[5] = (char *)&(test.dr)     - (char *)(&test);
   offsets[6] = (char *)&(test.wiph)   - (char *)(&test);
   offsets[7] = (char *)&(test.E_int_old) - (char *)(&test);
   offsets[8] = (char *)&(test.dV_old)  - (char *)(&test);
   offsets[9] = (char *)&(test.dE_cool) - (char *)(&test);

   MPI_Type_create_struct( count , blocksize , offsets , types , cell_mpi );
   MPI_Type_commit( cell_mpi );

}

void copy_cell_to_lite( struct cell * c , struct cell_lite * cl ){
  
   memcpy( cl->prim   , c->prim   , NUM_Q*sizeof(double) ); 
   memcpy( cl->cons   , c->cons   , NUM_Q*sizeof(double) ); 
   memcpy( cl->RKcons , c->RKcons , NUM_Q*sizeof(double) );
   memcpy( cl->grad   , c->grad   , NUM_Q*sizeof(double) );

   cl->riph   = c->riph;
   cl->dr     = c->dr;
   cl->wiph   = c->wiph; 

   cl->E_int_old   = c->E_int_old;
   cl->dV_old      = c->dV_old;
   cl->dE_cool     = c->dE_cool; 
 
}

void copy_lite_to_cell( struct cell_lite * cl , struct cell * c ){ 

   memcpy( c->prim   , cl->prim   , NUM_Q*sizeof(double) ); 
   memcpy( c->cons   , cl->cons   , NUM_Q*sizeof(double) ); 
   memcpy( c->RKcons , cl->RKcons , NUM_Q*sizeof(double) );
   memcpy( c->grad   , cl->grad   , NUM_Q*sizeof(double) );

   c->riph   = cl->riph;
   c->dr     = cl->dr;
   c->wiph   = cl->wiph;

   c->E_int_old = cl->E_int_old;
   c->dV_old    = cl->dV_old;
   c->dE_cool   = cl->dE_cool;


}

void generate_sendbuffer( struct domain * theDomain , struct cell_lite * pl , struct cell_lite * pr , int mode ){

   struct cell * theCells = theDomain->theCells;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   int i;

   for( i=0 ; i<Ng ; ++i ){
      if( mode==1 ){
         struct cell * cL = theCells+Ng+i;
         struct cell * cR = theCells+Nr-2*Ng+i;
         copy_cell_to_lite( cL , pl+i );
         copy_cell_to_lite( cR , pr+i );
      }else if( mode==2 ){
         struct cell * cL = theCells+i;
         struct cell * cR = theCells+Nr-Ng+i;
         if( rank != 0      ) copy_lite_to_cell( pl+i , cL );
         if( rank != size-1 ) copy_lite_to_cell( pr+i , cR );
      }
   }

}

void exchange_data( struct domain * theDomain ){

   MPI_Datatype cell_mpi = {0}; 
   generate_mpi_cell( &cell_mpi );

   int rank = theDomain->rank;
   int size = theDomain->size;

   int left = rank-1;
   if( left == -1 ) left = size-1;
   int right = rank+1;
   if( right == size ) right = 0;

   int Ng = theDomain->Ng;

   int tag = 0;
   MPI_Status status;

////////////
//Send Cells
////////////

   struct cell_lite pl_send[Ng];
   struct cell_lite pr_send[Ng];
   struct cell_lite pl_recv[Ng];
   struct cell_lite pr_recv[Ng];

//Build up list of cells to send...
   generate_sendbuffer( theDomain , pl_send , pr_send , 1 );
//Send!
   MPI_Sendrecv( pl_send , Ng , cell_mpi ,  left , tag+2,
                 pr_recv , Ng , cell_mpi , right , tag+2, MPI_COMM_WORLD , &status);
   MPI_Sendrecv( pr_send , Ng , cell_mpi , right , tag+3,
                 pl_recv , Ng , cell_mpi ,  left , tag+3, MPI_COMM_WORLD , &status);
//Now take the list of cells and put them into the appropriate locations...
   generate_sendbuffer( theDomain , pl_recv , pr_recv , 2 );

   MPI_Type_free( &cell_mpi );
}
