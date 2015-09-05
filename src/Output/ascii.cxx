
#include <fstream>
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
        struct cell * c = &(theCells[i]);
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
    fprintf(oFile, "mass loss:              %s \n",
             theDomain->mass_loss.c_str());


    if ( theDomain->SNe.size() > 0 )
    {
        fprintf(oFile, "Number of SNe:          %lu \n", theDomain->SNe.size());

        char SNe_filename[256] = "";
        strcat(SNe_filename, theDomain->output_prefix.c_str());
        strcat(SNe_filename, "SNe.dat");
        FILE * SNe_oFile = fopen(SNe_filename,"w");
        fprintf(SNe_oFile, "# SNe times [s]     ");
        fprintf(SNe_oFile, " initial mass [g]    ");
        fprintf(SNe_oFile, " ejecta mass [g]    ");
        fprintf(SNe_oFile, " ejecta mass (Z) [g] ");
        fprintf(SNe_oFile, " wind mass [g] \n");
        for (auto SN : theDomain->SNe)
        {
            fprintf(SNe_oFile, "%e         ", SN.lifetime);
            fprintf(SNe_oFile, "%e         ", SN.mass);
            fprintf(SNe_oFile, "%e         ", SN.mass_ejecta);
            fprintf(SNe_oFile, "%e         ", SN.mass_ejecta_Z);
            fprintf(SNe_oFile, "%e         ", SN.mass_winds);
            fprintf(SNe_oFile, "\n");
        }

        fclose(SNe_oFile);
    }

    time_t current_time = time(NULL);
    fprintf(oFile, "Created at: %s \n", ctime(&current_time));

    fclose(oFile);

    return 0;
}

std::size_t count_lines_in_file( const std::string filename )
{

    if ( !fs::exists(filename) )
    {
        std::cerr << "Error: File (\"" << filename 
                  << "\" doesn't exist. Can't count lines." << std::endl;
        return 0;
    }

   std::size_t num_lines = 0;
   std::string line;
   std::ifstream in_File(filename);   
   while (std::getline(in_File , line))
        ++num_lines;

return num_lines;
}
