#include <iostream>
#include <assert.h>
#include <time.h>

extern "C" {
#include <grackle.h>
}

#include "ascii.H"
#include "../structure.H"
#include "../geometry.H"

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;


void output( struct domain * theDomain , const char * filestart , 
             const double t )
{

    struct cell * theCells = theDomain->theCells;
    int Nr = theDomain->Nr;

    char filename[256] = "";
    strcat(filename, theDomain->output_prefix.c_str());
    strcat(filename, filestart);
    strcat(filename, ".dat");

    FILE * pFile = fopen( filename , "w" );
    fprintf(pFile,"# time = %le [s] \n", t);
    // fprintf(pFile,"# r           dr           dV           Density      Pressure     Velocity     Z\n");
    fprintf(pFile,"# r                  dr                 dV                 Density            Pressure           Velocity           Z\n");

    int i_min = 0;
    int i_max = Nr;

    int i,q;
    for( i=i_min ; i<i_max ; ++i )
    {
        struct cell * c = theCells+i;
        double rp = c->riph;
        double dr = c->dr; 
        double rm = rp-dr;
        double dV = get_dV( rp , rm );
        fprintf(pFile,"%18.10e %18.10e %18.10e ",rp,dr,dV);
        for( q=0 ; q<NUM_Q ; ++q )
        {
            fprintf(pFile,"%18.10e ",c->prim[q]);
        }
        #ifndef NDEBUG
        if(c->prim[PPP] < theDomain->theParList.Pressure_Floor)
        {
            printf("-------ERROR--------- \n");
            printf("In output() \n");
            printf("c->prim[PPP] = %e at i=%d \n", c->prim[PPP], i);
            printf("Pressure floor should be = %e \n", theDomain->theParList.Pressure_Floor);
            assert(0);
        }
        #endif
        fprintf(pFile,"\n");
    }
    fclose( pFile );

}

int overview( struct domain * theDomain )
{
    // prints an overview of key parameters into a datafile


    // some of this should be moved into a seperate file
    // which changes depending on the initial conditions / which run we're doing

    char overview_filename[256] = "";
    strcat(overview_filename, theDomain->output_prefix.c_str());
    strcat(overview_filename, "overview.dat");

    fs::path overview_path(overview_filename);
    if ( fs::exists(overview_path) )
    {
        std::cerr << "Warning: overview already exists; not overwriting" 
                  << std::endl;
        if ( theDomain->theParList.ICs.compare("restart") == 0 )
        {
            // if we're just restarting, we expect overview to already exist
            return 0;
        }
        else
        {
            // the overview file already existing might lead to some problems
            // exit now, and figure out what's happening.
            return 1;
        }
    }

    FILE * oFile = fopen(overview_filename,"w");

    fprintf(oFile, "Metallicity:            %e \n",
             theDomain->metallicity);
    fprintf(oFile, "Background Density:     %e \n",
             theDomain->background_density);
    fprintf(oFile, "Background Temperature: %e \n",
             theDomain->background_temperature);
    fprintf(oFile, "With cooling:           %d \n",
             theDomain->theParList.With_Cooling);
    fprintf(oFile, "Cluster Mass [M_sol]:   %e \n",
             theDomain->cluster_mass);
    fprintf(oFile, "seed:                   %d \n",
             theDomain->seed);

    if ( theDomain->SNe_times.size() > 0 )
    {
        fprintf(oFile, "Number of SNe:       %lu \n", theDomain->SNe_times.size());

        char SNe_times_filename[256] = "";
        strcat(SNe_times_filename, theDomain->output_prefix.c_str());
        strcat(SNe_times_filename, "SNe_times.dat");
        FILE * SNe_times_oFile = fopen(SNe_times_filename,"w");
        fprintf(SNe_times_oFile, "# SNe times [s]: \n");
        for (auto SNe_time : theDomain->SNe_times)
            fprintf(SNe_times_oFile, "%e \n", SNe_time);
        fclose(SNe_times_oFile);
    }

    time_t current_time = time(NULL);
    fprintf(oFile, "Created at: %s \n", ctime(&current_time));

    fclose(oFile);

    return 0;
}
