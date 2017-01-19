
#include <iostream>
#include <string>
#include <stdexcept>
#include "../structure.H"

#include "initial_conditions.H"
#include "restart_ICs.H"
#include "../misc.H" // calc_dr
#include "../boundary.H"
#include "../Hydro/euler.H" // prim2cons, cons2prim, mindt
#include "../geometry.H"
#include "../blast.H"
#include "../Output/ascii.H" // count_lines_in_file, read_SNe


#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

const std::string Restart_ICs::class_name = "restart";


Restart_ICs::Restart_ICs() : Initial_Conditions( class_name )
{
    NL  = 0;
    rr  = nullptr;
    rho = nullptr;
    Pp  = nullptr;
    vr  = nullptr;
    Z   = nullptr;

    restart_filename = std::string("");
    restart_id = std::string("");

    time_restart = 0.0;
    last_finished_checkpoint = 0;
}

bool Restart_ICs::trust_LogZoning_flag() const
{
    return false;
}

int Restart_ICs::setICparams( struct domain * theDomain ,
                              const Mass_Loss * mass_loss )
{
    NL = this->get_table( restart_filename ); 
    if ( NL < 0 )
    {
        std::cerr << "Error reading in restart file" << std::endl;
        return 1;
    }


    std::cout << std::endl;
    std::cout << "Restart values: " << std::endl;


    // Overwrite using command line args
    int N_chk_tmp = theDomain->N_chk;
    if ( N_chk > 0 )
    {
        N_chk_tmp = N_chk;
        std::cout << "Overwriting number of checkpoints to: " << N_chk << std::endl;
    }

    double delta_time_tmp = theDomain->t_fin - theDomain->t_init;
    if ( delta_time > 0 )
    {
        delta_time_tmp = delta_time;
        std::cout << "Overwriting delta_time to: " << delta_time << std::endl;
    }

    if ( CFL > 0 )
    {
        theDomain->theParList.CFL = CFL;
        std::cout << "Overwriting CFL to: " << CFL << std::endl;
    }

    if (Cooling_Redshift > 0 )
    {
        theDomain->theParList.Cooling_Redshift = Cooling_Redshift;
        std::cout << "Overwriting Cooling_Redshift to: " << Cooling_Redshift << std::endl;
    }

    // Set values

    theDomain->t      = time_restart;
    theDomain->t_init = time_restart;
    theDomain->t_fin  = time_restart + delta_time_tmp;
    theDomain->nchk_0 = last_finished_checkpoint ;
    theDomain->nchk   = last_finished_checkpoint ;
    theDomain->N_chk  = N_chk_tmp + 1; // since we don't do a checkpoint when starting

    this->set_output_prefix( theDomain );
    this->add_SNe( theDomain , mass_loss );
    this->set_times( theDomain );

    std::cout << "time restart:   " << time_restart << std::endl;
    std::cout << "last_finished_checkpoint: " << last_finished_checkpoint << std::endl;
    std::cout << "restarting uuid: " << theDomain->output_prefix << std::endl;
    std::cout << std::endl;

    return 0;
}

void Restart_ICs::add_SNe( struct domain * theDomain ,
                           const Mass_Loss * mass_loss )
{
    std::string SNe_filename (restart_filename);

    std::vector<supernova> SNe = read_SNe(theDomain->output_prefix + "SNe.dat");
    
    assert( std::is_sorted( SNe.rbegin(),
                            SNe.rend(),
                            sort_by_lifetime) );

    const double tolerance = 1e-5;
    while ( (SNe.size() > 0) && 
            (SNe.back().lifetime <= theDomain->t*(1+tolerance)) )
    {
        SNe.pop_back();
    }

    theDomain->SNe = SNe;
}

void Restart_ICs::initial( double * prim , double r )
{
    // since we don't want any hard-coded initial conditions,
    // we don't want to use this method.

    // instead, everything should be set within setup_grid

    return;
}

int Restart_ICs::parse_command_line_args (  struct domain * theDomain , 
                                            int argc , 
                                            char * argv [] )
{
    std::string partial_ID("");
    if ( argc > 2 )
    {
        partial_ID = std::string(argv[2]);
    }  

    restart_filename = this->get_restart_filename(partial_ID);

    N_chk = -1;
    if ( argc > 3 )
    {
        N_chk = std::stoi(argv[3]);
    }  

    delta_time = -1;
    if ( argc > 4 )
    {
        delta_time = std::stod(argv[4]);
    }  

    CFL = -1;
    if ( argc > 5 )
    {
        CFL = std::stod(argv[5]);
    }  

    Cooling_Redshift = -1;
    if ( argc > 6 )
    {
        Cooling_Redshift = std::stod(argv[6]);
    }  

    return 0;
}




std::string Restart_ICs::get_restart_filename( const std::string partial_ID )
{
    // there's probably a better way to do this.
    // boost::regex perhaps?

    int highest_checkpoint_seen = -1;

    fs::path current_dir("."); // make this an input parameter?

    fs::directory_iterator end_itr;

    std::string output_prefix;

    for ( fs::directory_iterator dir_itr{current_dir} ;
          dir_itr != end_itr ; 
          ++dir_itr )
    {
        std::string tmp_filename = dir_itr->path().string();

        if (tmp_filename.find(partial_ID) == std::string::npos)
        {
            continue;
        }
        std::size_t pos_left = tmp_filename.find("checkpoint_");
        if (tmp_filename.find("checkpoint_") == std::string::npos )
        {
            continue;
        }
        pos_left += std::string("checkpoint_").size();

        std::size_t pos_right = tmp_filename.find(".dat");

        int checkpoint = std::stoi(tmp_filename.substr(pos_left, 
                                                       pos_right-pos_left));

        std::string tmp_output_prefix = this->filename_to_prefix(tmp_filename);
        if ( output_prefix.empty() )
        {
            output_prefix = tmp_output_prefix;
        }

        if ( output_prefix.compare(tmp_output_prefix) != 0)
        {
            throw std::runtime_error(std::string("partial_ID: '") 
                                     + partial_ID
                                     + std::string("' not sufficiently unique"));
        }

        if ( checkpoint > highest_checkpoint_seen )
        {
            restart_filename.swap(tmp_filename);
            highest_checkpoint_seen = checkpoint;
        }

    }

    // check that file exists
    fs::path restart_path(restart_filename);
    assert(fs::exists(restart_path));

    if ( highest_checkpoint_seen >= 0 )
    {
        last_finished_checkpoint = highest_checkpoint_seen;
    }

    return restart_filename;

}


int Restart_ICs::get_table( const std::string filename )
{
    int nL = count_lines_in_file(filename) - 2;
    if ( nL < 0 ) return nL;
    rr  = (double *) malloc( nL*sizeof(double) );
    double dr;
    double dV;
    rho = (double *) malloc( nL*sizeof(double) );
    Pp  = (double *) malloc( nL*sizeof(double) );
    vr  = (double *) malloc( nL*sizeof(double) );
    Z   = (double *) malloc( nL*sizeof(double) );
    FILE * pFile = fopen(filename.c_str(),"r");
    char tmp[1024];
    fscanf(pFile,"%s %s %s %le %s \n",
    tmp, tmp, tmp, &time_restart, tmp);

    fgets(tmp, sizeof(tmp), pFile); // header line

    for( int l=0 ; l<nL ; ++l )
    {
        fscanf(pFile,"%le %le %le %le %le %le %le\n",
                &(rr[l]),&dr,&dV,&(rho[l]),&(Pp[l]),&(vr[l]),&(Z[l]));
    }
    fclose(pFile);
    return nL;
}



void Restart_ICs::setup_grid( struct domain * theDomain )
{
    // this is a combination of the setup_grid and setup_cells,
    // adapted from the base class Initial_Conditions
    theDomain->Ng = NUM_G;

    const int Nr = NL;
    theDomain->Nr = Nr;

    theDomain->theCells = (struct cell *) malloc( Nr*sizeof(struct cell));

    for( int i=0 ; i<Nr ; ++i )
    {
        struct cell * c = &(theDomain->theCells[i]);
        c->riph         = rr[ i];
        c->prim[RHO]    = rho[i];
        c->prim[PPP]    = Pp[ i];
        c->prim[VRR]    = vr[ i];
        c->prim[ZZZ]    = Z[  i];
    }
    calc_dr( theDomain );
    set_wcell( theDomain );

    for( int i=0 ; i<Nr ; ++i )
    {
        struct cell * c = &(theDomain->theCells[i]);
        double rp = c->riph;
        double rm = rp - c->dr;
        double dV = get_dV( rp , rm );
        prim2cons( c->prim , c->cons , dV );
        cons2prim( c->cons , c->prim , dV );
        c->E_int_old = E_int_from_cons( c->cons );
        c->dV_old = dV;
    }

    boundary( theDomain );
}

void Restart_ICs::set_times( struct domain * theDomain )
{
    // just use the times which were set in setupDomain
    // (i.e. just evolve for the time given in the parameter file)

    return;

}

void Restart_ICs::set_output_prefix( struct domain * theDomain )
{
    theDomain->output_prefix = this->filename_to_prefix(restart_filename);
}

std::string Restart_ICs::filename_to_prefix( const std::string filename ) const
{

    if (filename.find("checkpoint_") != std::string::npos )
    {
        const std::string basename = fs::basename(filename);
        const unsigned int id_start = basename.find("checkpoint_");

        return basename.substr(0,id_start);
    }
    else
    {
        throw std::runtime_error(std::string("Couldn't find an output prefix using filename: ")
                                 + filename);
    }

}


void Restart_ICs::free_table()
{
    free(rr);
    free(rho);
    free(Pp);
    free(vr);
    free(Z);
}

Restart_ICs::~Restart_ICs()
{
    this->free_table();
}

