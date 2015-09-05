
#include <iostream>
#include <string>
#include <cmath>

#include <uuid/uuid.h>
extern "C" {
#include <grackle.h>
}

#include <assert.h>

#include "blast.H"
#include "constants.H"
#include "cooling.H"
#include "domain.H"
#include "geometry.H"
#include "riemann.H"
#include "structure.H"
#include "Initial/initial_conditions.H"
#include "Hydro/euler.H"
#include "Output/ascii.H"

#include <algorithm>


int setupDomain( struct domain * theDomain , 
                 Initial_Conditions * ICs , 
                 Mass_Loss * mass_loss )
{

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

    theDomain->nrpt   = -1;
    theDomain->nsnp   = -1;
    theDomain->nchk   = -1;
    theDomain->nchk_0 =  0;

    uuid_t  id_binary;
    char id_ascii[36];
    uuid_generate(id_binary);
    uuid_unparse(id_binary, id_ascii);
    printf("generated uuid: %s \n", id_ascii);
    theDomain->output_prefix = std::string(id_ascii).append("_");

    error = ICs->setICparams( theDomain );
    if ( error==1 ) return(error);
    setHydroParams( theDomain );
    setRiemannParams( theDomain );

    // everything below should probably go in a separate function?

    assert( std::is_sorted( theDomain->SNe.rbegin(),
                            theDomain->SNe.rend(),
                            sort_by_lifetime) );

    ICs->set_times( theDomain );

    return 0;

}


void freeDomain( struct domain * theDomain ){
    free( theDomain->theCells );
}


void check_dt( struct domain * theDomain , double * dt ){

    double t = theDomain->t;
    double tmax = theDomain->t_fin;
    int final=0;
    if( t + *dt > tmax )
    {
        *dt = tmax-t;
        final=1;
    }

    FILE * abort = NULL;
    abort = fopen("abort","r");
    if( abort ){ final = 1; fclose(abort); }

    if( final ) theDomain->final_step = 1;

}

// void snapshot( struct domain * , char * );

void possiblyOutput( struct domain * theDomain , int override )
{

    double t     = theDomain->t;
    double t_min = theDomain->t_init;
    double t_fin = theDomain->t_fin;
    int Nrpt     = theDomain->N_rpt;
    // int Nsnp = theDomain->N_snp;
    int Nchk     = theDomain->N_chk;
    int nchk_0   = theDomain->nchk_0;

    int LogOut = theDomain->theParList.Out_LogTime;
    int n0;

    n0 = static_cast<int> ( (t-t_min)*Nrpt / (t_fin-t_min) ) ;
    if( LogOut ) n0 = static_cast <int> ( Nrpt*std::log(t/t_min)
                                              /std::log(t_fin/t_min) );
    n0 += nchk_0;
    if( theDomain->nrpt < n0 || override )
    {
        theDomain->nrpt = n0;
        printf("t = %.3e\n",t);
    }

    n0 = static_cast<int> ( (t-t_min)*Nchk / (t_fin-t_min) );
    if( LogOut ) n0 = static_cast<int> ( Nchk*std::log(t/t_min)
                                             /std::log(t_fin/t_min) );
    n0 += nchk_0;
    if( (theDomain->nchk < n0 && Nchk>0) || override )
    {
        theDomain->nchk = n0;
        char filename[256];
        if( !override )
        {
            printf("Creating Checkpoint #%04d...\n",n0);
            sprintf(filename,"checkpoint_%04d",n0);
            output( theDomain , filename, t );
        }
        else
        {
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


