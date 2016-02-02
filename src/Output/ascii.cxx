
#include <fstream>
#include <iostream>
#include <assert.h>
#include <time.h>
#include <stdexcept>
#include <sstream>

extern "C" {
#include <grackle.h>
}

#include "ascii.H"
#include "../structure.H"
#include "../geometry.H"
#include "../blast.H" // sort_by_lifetime

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;


void create_checkpoint( struct domain * theDomain , const char * filestart , 
                 const double t )
{

    const struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;

    char filename[256] = "";
    strcat(filename, theDomain->output_prefix.c_str());
    strcat(filename, filestart);
    strcat(filename, ".dat");

    FILE * pFile = fopen( filename , "w" );
    fprintf(pFile,"# time = %le [s] \n", t);
    fprintf(pFile,"# r                  dr                 dV                 Density            Pressure           Velocity           Z\n");

    const int i_min = 0;
    const int i_max = Nr;

    for( int i=i_min ; i<i_max ; ++i )
    {
        const struct cell * c = &(theCells[i]);
        const double rp = c->riph;
        const double dr = c->dr; 
        const double rm = rp-dr;
        const double dV = get_dV( rp , rm );
        fprintf(pFile,"%18.10e %18.10e %18.10e ",rp,dr,dV);
        for( int q=0 ; q<NUM_Q ; ++q )
        {
            fprintf(pFile,"%18.10e ",c->prim[q]);
        }
        #ifndef NDEBUG
        if(c->prim[PPP] < theDomain->theParList.Pressure_Floor)
        {
            printf("-------ERROR--------- \n");
            printf("In create_checkpoint() \n");
            printf("c->prim[PPP] = %e at i=%d \n", c->prim[PPP], i);
            printf("Pressure floor should be = %e \n", theDomain->theParList.Pressure_Floor);
            assert(0);
        }
        #endif
        fprintf(pFile,"\n");
    }
    fclose( pFile );

}

void overviews(  struct domain * theDomain ,
                Mass_Loss * mass_loss ,
                Cooling * cooling )
{

    main_overview( theDomain , mass_loss, cooling );

    SNe_overview( theDomain );

    inputs_overview( theDomain , mass_loss, cooling );


}

void main_overview(  struct domain * theDomain ,
                    Mass_Loss * mass_loss ,
                    Cooling * cooling )
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
            return;
        }
        else
        {
            // the overview file already existing might lead to some problems
            // exit now, and figure out what's happening.
            throw std::runtime_error("Overview already exists");
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
             cooling->with_cooling);
    fprintf(oFile, "Cooling Type:           %s \n",
             cooling->name.c_str() );
    fprintf(oFile, "Cluster Mass [M_sol]:   %e \n",
             theDomain->cluster_mass);
    fprintf(oFile, "seed:                   %d \n",
             theDomain->seed);
    fprintf(oFile, "mass loss:              %s \n",
            mass_loss->name.c_str());

    fprintf(oFile, "Number of SNe:          %lu \n", theDomain->SNe.size());

    time_t current_time = time(NULL);
    fprintf(oFile, "Created at: %s \n", ctime(&current_time));

    fclose(oFile);

    return;
}

void SNe_overview( struct domain * theDomain )
{
    // prints the initial state of the SNe vector into a data file

    char SNe_filename[256] = "";
    strcat(SNe_filename, theDomain->output_prefix.c_str());
    strcat(SNe_filename, "SNe.dat");

    fs::path SNe_path(SNe_filename);
    if ( fs::exists(SNe_path) )
    {
        std::cerr << "Warning: SNe overview already exists; not overwriting" 
                  << std::endl;
        if ( theDomain->theParList.ICs.compare("restart") == 0 )
        {
            // if we're just restarting, we expect SNe overview to already exist
            return;
        }
        else
        {
            // the SNe overview file already existing might lead to some problems
            // exit now, and figure out what's happening.
            throw std::runtime_error("SNe Overview already exists");
        }
    }

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

    return;
}

void inputs_overview( struct domain * theDomain ,
                      Mass_Loss * mass_loss ,
                      Cooling * cooling )
{

    // prints a copy of the input parameters into a data file


    // some of this should be moved into a seperate file
    // which changes depending on the initial conditions / which run we're doing

    char inputs_overview_filename[256] = "";
    strcat(inputs_overview_filename, theDomain->output_prefix.c_str());
    strcat(inputs_overview_filename, "inputs.dat");

    fs::path inputs_overview_path(inputs_overview_filename);
    if ( fs::exists(inputs_overview_path) )
    {
        std::cerr << "Warning: inputs overview already exists; not overwriting" 
                  << std::endl;
        if ( theDomain->theParList.ICs.compare("restart") == 0 )
        {
            // if we're just restarting, we expect overview to already exist
            return;
        }
        else
        {
            // the overview file already existing might lead to some problems
            // exit now, and figure out what's happening.
            throw std::runtime_error("Input Overview already exists");
        }
    }

    FILE * oFile = fopen(inputs_overview_filename,"w");


    fprintf(oFile, "T_Start:          %e \n",
             theDomain->theParList.t_min);
    fprintf(oFile, "T_End:            %e \n",
             theDomain->theParList.t_max);


    fprintf(oFile, "Num_Reports:      %d \n",
             theDomain->theParList.NumRepts);
    fprintf(oFile, "Num_Checkpoints:  %d \n",
             theDomain->theParList.NumChecks);
    fprintf(oFile, "Use_Logtime:      %d \n",
             theDomain->theParList.Out_LogTime);

    fprintf(oFile, "\n");


    fprintf(oFile, "Num_R:            %d \n",
             theDomain->theParList.Num_R);
    fprintf(oFile, "R_Min:            %e \n",
             theDomain->theParList.rmin);
    fprintf(oFile, "R_Max:            %e \n",
             theDomain->theParList.rmax);
    fprintf(oFile, "Log_Zoning:       %d \n",
             theDomain->theParList.LogZoning);
    fprintf(oFile, "Log_Radius:       %e \n",
             theDomain->theParList.LogRadius);

    fprintf(oFile, "\n");


    fprintf(oFile, "CFL:              %e \n",
             theDomain->theParList.CFL);
    fprintf(oFile, "PLM:              %d \n",
             theDomain->theParList.PLM);
    fprintf(oFile, "RK2:              %d \n",
             theDomain->theParList.RK2);
    fprintf(oFile, "H_0:              %e \n",
             theDomain->theParList.H_0);
    fprintf(oFile, "H_1:              %e \n",
             theDomain->theParList.H_1);
    fprintf(oFile, "Riemann_Solver:   %d \n", 
             theDomain->theParList.Riemann_Solver);
    fprintf(oFile, "Mesh_Motion:      %d \n",
             theDomain->theParList.Mesh_Motion);
    fprintf(oFile, "Density_Floor:    %e \n",
             theDomain->theParList.Density_Floor);
    fprintf(oFile, "Pressure_Floor:   %e \n",
             theDomain->theParList.Pressure_Floor);

    fprintf(oFile, "\n");


    fprintf(oFile, "With_Cooling:     %d \n",
             theDomain->theParList.With_Cooling);
    fprintf(oFile, "Cooling_Type:     %s \n",
             cooling->name.c_str());

    fprintf(oFile, "\n");

    fprintf(oFile, "Adiabatic_Index:  %e \n",
             theDomain->theParList.Adiabatic_Index);

    fprintf(oFile, "\n");

    fprintf(oFile, "ICs:              %s \n",
             theDomain->theParList.ICs.c_str());

    fprintf(oFile, "\n");

    fprintf(oFile, "mass_loss:        %s \n",
            mass_loss->name.c_str());

    fclose(oFile);


    return;
}


std::size_t count_lines_in_file( const std::string filename )
{

    if ( !fs::exists(filename) )
    {
        std::cerr << "Error: File \"" << filename 
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

std::vector<supernova> read_SNe( const std::string filename)
{

    std::vector<supernova> SNe;
    supernova SN_tmp;

    const int nL = count_lines_in_file(filename) - 1;
    if ( nL < 0 ) return SNe;

    double SN_time;
    double SN_mass;
    double SN_mass_ejecta;
    double SN_mass_ejecta_Z;
    double SN_mass_winds;


    FILE * pFile = fopen(filename.c_str(),"r");
    char tmp[1024];
    fgets(tmp, sizeof(tmp), pFile); // header line

    for( int l=0 ; l<nL ; ++l )
    {
        fscanf(pFile,"%le %le %le %le %le\n",
                &SN_time, &SN_mass, 
                &SN_mass_ejecta, &SN_mass_ejecta_Z,
                &SN_mass_winds );

        SN_tmp.mass             = SN_mass;
        SN_tmp.mass_ejecta      = SN_mass_ejecta;
        SN_tmp.mass_ejecta_Z    = SN_mass_ejecta_Z;
        SN_tmp.mass_winds       = SN_mass_winds;
        SN_tmp.lifetime         = SN_time;

        SNe.push_back(SN_tmp);
    }
    fclose(pFile);

    std::sort(SNe.rbegin(), SNe.rend(), sort_by_lifetime);

    return SNe;
}



