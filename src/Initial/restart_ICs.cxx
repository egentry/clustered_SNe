
#include <iostream>
#include <string>
#include "../structure.H"

#include "initial_conditions.H"
#include "restart_ICs.H"
#include "../misc.H" // calc_dr
#include "../boundary.H"
#include "../Hydro/euler.H" // prim2cons, cons2prim, mindt
#include "../geometry.H"
#include "../blast.H"
#include "../Output/ascii.H"


#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

Restart_ICs::Restart_ICs()
{
    NL  = 0;
    rr  = NULL;
    rho = NULL;
    Pp  = NULL;
    vr  = NULL;
    Z   = NULL;

    restart_filename = std::string("init");
    restart_id = std::string("");

    time_restart = 0.0;
    checkpoints_finished = 0;
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

    const double delta_time = theDomain->t_fin - theDomain->t_init;

    theDomain->t      = time_restart;
    theDomain->t_init = time_restart;
    theDomain->t_fin  = time_restart + delta_time;
    theDomain->nchk_0 = checkpoints_finished ;
    theDomain->nchk   = checkpoints_finished ;
    theDomain->N_chk += 1; // since we don't do a checkpoint when starting

    this->set_output_prefix( theDomain );
    this->add_SNe( theDomain , mass_loss );
    this->set_times( theDomain );

    std::cout << std::endl;
    std::cout << "Restart values:" << std::endl;
    std::cout << "time restart: " << time_restart << std::endl;
    std::cout << "checkpoints_finished: " << checkpoints_finished << std::endl;
    std::cout << "output_prefix: " << theDomain->output_prefix << std::endl;
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

    
    while ( (SNe.size() > 0) && 
            (SNe.back().lifetime < theDomain->t) )
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

    return 0;
}




std::string Restart_ICs::get_restart_filename( const std::string partial_ID )
{
    // there's probably a better way to do this.
    // boost::regex perhaps?

    int highest_checkpoint_seen = -1;

    fs::path current_dir("."); // make this an input parameter?

    fs::directory_iterator end_itr;

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
        checkpoints_finished = highest_checkpoint_seen;
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
        // if(l==1)
        // {
        //    printf("r = %le \n", rr[l]);
        //    printf("dr = %le \n", dr);
        //    printf("dV = %le \n", dV);
        //    printf("rho = %le \n", rho[l]);
        //    printf("P = %le \n", Pp[l]);
        //    printf("vr = %le \n", vr[l]);
        //    printf("Z = %le \n", Z[l]);         
        // }

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
        c->E_int_old = E_int_from_cons( c );
        c->dV_old = dV;
    }

    boundary( theDomain );
}

void Restart_ICs::set_times( struct domain * theDomain )
{
    // do basically the same thing as the base class,
    // except don't change t or t_init

    assert( std::is_sorted( theDomain->SNe.rbegin(),
                            theDomain->SNe.rend(),
                            sort_by_lifetime) );

    if ( theDomain->SNe.size() > 0 )
    {

        double t_last_SN  = theDomain->SNe.front().lifetime;

        theDomain->t_fin  += t_last_SN;

    }

}

void Restart_ICs::set_output_prefix( struct domain * theDomain )
{

    if (restart_filename.find("checkpoint_") != std::string::npos )
    {
        const std::string restart_basename = fs::basename(restart_filename);
        const unsigned int id_start = restart_basename.find("checkpoint_");

        theDomain->output_prefix = restart_basename.substr(0,id_start);
    }
    else
    {
        std::cerr << "Couldn't find an output prefix using restart_filename: "
                  << restart_filename << std::endl;
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


