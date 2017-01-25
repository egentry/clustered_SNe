
#include <iostream>
#include <string>
#include <cmath> // std::abs
#include <uuid/uuid.h>
#include <stdexcept>



#include "../Hydro/euler.H" // prim2cons, cons2prim
#include "../boundary.H"
#include "../geometry.H" // get_moment_arm, get_dV
#include "../misc.H" // calc_dr, E_int_from_*
#include "../blast.H"

#include "initial_conditions.H"
#include "chevalier_ICs.H"
#include "cluster_SNe_ICs.H"
#include "ejecta_ICs.H"
#include "isentropic_ICs.H"
#include "messy_ICs.H"
#include "restart_ICs.H"
#include "shocktube_ICs.H"
#include "Thornton_parameter_study_ICs.H"
#include "uniform_ICs.H"
#include "conduction_study_ICs.H"
#include "Sedov_ICs.H"

#include <boost/algorithm/string.hpp>    

Initial_Conditions * select_initial_conditions( std::string IC_name )
{

    std::cout << "Using ICs: " << IC_name << std::endl;

    boost::algorithm::to_lower(IC_name);


    if ( IC_name.compare(boost::to_lower_copy(Chevalier_ICs::class_name)) == 0 )
    {
        return new Chevalier_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Ejecta_ICs::class_name)) == 0 )
    {
        return new Ejecta_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Cluster_SNe_ICs::class_name)) == 0 )
    {
        return new Cluster_SNe_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Isentropic_ICs::class_name)) == 0 )
    {
        return new Isentropic_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Messy_ICs::class_name)) == 0 )
    {
        return new Messy_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Restart_ICs::class_name)) == 0 )
    {
        return new Restart_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Shocktube_ICs::class_name)) == 0 )
    {
        return new Shocktube_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Thornton_Parameter_Study_ICs::class_name)) == 0 )
    {
        return new Thornton_Parameter_Study_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Uniform_ICs::class_name)) == 0 )
    {
        return new Uniform_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Conduction_Study_ICs::class_name)) == 0 )
    {
        return new Conduction_Study_ICs;
    }

    if ( IC_name.compare(boost::to_lower_copy(Sedov_ICs::class_name)) == 0 )
    {
        return new Sedov_ICs;
    }

    throw std::invalid_argument("Initial condition name didn't match known names");

}

//************ Initial_Conditions (abstract base class) *********************//

Initial_Conditions::Initial_Conditions( const std::string name ) : name(name)
{}

int Initial_Conditions::parse_command_line_args( struct domain * theDomain , 
                                                 int argc , 
                                                 char * argv [] ){ return 0; }

bool Initial_Conditions::trust_LogZoning_flag() const
{
    return true;
}

void Initial_Conditions::add_SNe( struct domain * theDomain ,
                                  const Mass_Loss * mass_loss )
{
    return;
}

int Initial_Conditions::setICparams( struct domain * theDomain ,
                                     const Mass_Loss * mass_loss)
{
    this->set_output_prefix( theDomain );
    this->set_times( theDomain );

    return 0;

}

void Initial_Conditions::setup_cells( struct domain * theDomain )
{

    int i;
    struct cell * theCells = theDomain->theCells;
    int Nr = theDomain->Nr;

    for( i=0 ; i<Nr ; ++i )
    {
        struct cell * c = &(theCells[i]);
        double rp = c->riph;
        double rm = rp - c->dr;
        c->wiph = 0.0; 
        double r = get_moment_arm( rp , rm );
        this->initial( c->prim , r ); 
        double dV = get_dV( rp , rm );
        prim2cons( c->prim , c->cons , dV );
        cons2prim( c->cons , c->prim , dV );
        c->E_int_old = E_int_from_cons( c->cons );
        c->dV_old = dV;

        c->multiphase = 0;
        for ( int q=0 ; q<NUM_Q ; ++q)
        {
            c->prim_hot[q]  = 0;
            c->prim_cold[q] = 0;
            c->cons_hot[q]  = 0;
            c->cons_cold[q] = 0;

            c->RKcons_hot[q]=0;
            c->RKcons_cold[q]=0;
        }
        c->V_hot  = 0;
        c->V_cold = 0;

        c->x_hot  = 0;
        c->x_cold = 0;
        c->y_hot  = 0;
        c->y_cold = 0;
        c->z_hot  = 0;
        c->z_cold = 0;

        c->E_kin_initial = 0;
        c->E_int_initial = 0;
    }

    boundary( theDomain );

}

void Initial_Conditions::setup_grid( struct domain * theDomain )
{
    theDomain->Ng = NUM_G;
    int Num_R = theDomain->theParList.Num_R;
    int LogZoning = theDomain->theParList.LogZoning;

    double Rmin = theDomain->theParList.rmin;
    double Rmax = theDomain->theParList.rmax;

    int Nr = Num_R;

    theDomain->Nr = Nr;
    theDomain->theCells = (struct cell *) malloc( Nr*sizeof(struct cell));

    int i;

    double dx = 1.0/static_cast<double>(Num_R);
    double x0 = 0;
    double R0 = theDomain->theParList.LogRadius;
    for( i=0 ; i<Nr ; ++i )
    {
        double xp = x0 + (i+1.0)*dx;
        double rp;
        if( LogZoning == 0 )
        {
            rp = Rmin + xp*(Rmax-Rmin);
        }
        else if( LogZoning == 1 )
        {
            rp = Rmin*pow(Rmax/Rmin,xp);
        }
        else
        {
            rp = R0*pow(Rmax/R0,xp) + Rmin-R0 + (R0-Rmin)*xp;
        }
        theDomain->theCells[i].riph = rp;
    }
    calc_dr( theDomain );

    this->setup_cells( theDomain );
};

void Initial_Conditions::possibly_extend_grid( struct domain * theDomain )
{

    // ============================================= //
    //
    //  Checks if more cells need to be added to the grid
    //
    //  Inputs:
    //     - theDomain - the standard domain struct used throughout
    //     - extend_fraction - Optional, default = 0.25
    //                       - If we currently have N cell,s
    //                       - this will create and additional:
    //                         N * extend_fraction of cells;
    //
    //  Returns:
    //    void
    //
    //  Side effects:
    //     - Inherits side effects of extend_grid
    //
    //  Notes:
    //    - Assumes shock is travelling outwards (extends outer boundary)
    //    - Uses shock-finding to determine if grid should be extended
    //    - This will throw off energy conservation calculations in analysis
    // ============================================= //

    const unsigned int Nr = theDomain->Nr;
    const double R_max    = theDomain->theCells[Nr-2].riph;
    const double R_shock  = theDomain->R_shock;


    const double radius_threshold = .9;
    if ( (R_shock / R_max) > radius_threshold )
    {
        this->extend_grid( theDomain );
    }
}

double Initial_Conditions::find_shock( const struct domain * theDomain ) const
{
    // ============================================= //
    //
    //  Finds the location of the shock
    //
    //  Inputs:
    //     - theDomain - the standard domain struct used throughout
    //
    //  Returns:
    //    - R_shock - outer boundary of first unshocked zone
    //
    //  Side effects:
    //     - None
    //
    //  Notes:
    //    - Assumes shock is travelling outwards (extends outer boundary)
    //    - Assumes constant density background
    //    - Used for determining when to:
    //          - Extend grid
    //          - Use previously cached cooling results
    // ============================================= //

    const unsigned int Nr = theDomain->Nr;
    const double background_density = theDomain->theCells[Nr-2].prim[RHO];

    const double density_tolerance = 1e-2;

    int i_shock = 0;
    double R_shock = 0;
    for (int i = 1 ; i < Nr ; ++i )
    {
        const struct cell * c = &(theDomain->theCells[i]);
        if ( ((c->prim[RHO] - background_density) / background_density) 
                > density_tolerance )
        {
            i_shock = i;
            R_shock = c->riph;
        }
    }

    return R_shock;
}


void Initial_Conditions::extend_grid( struct domain * theDomain, 
                                        const double extend_fraction )
{
    // ============================================= //
    //
    //  Adds an additional 'extend_fraction' more cells 
    //  to the outside of the grid
    //
    //  Inputs:
    //     - theDomain - the standard domain struct used throughout
    //     - extend_fraction - Optional, default = 0.25
    //                       - If we currently have N cell,s
    //                       - this will create and additional:
    //                         N * extend_fraction of cells;
    //
    //  Returns:
    //    void
    //
    //  Side effects:
    //     - changes domain size; adds new cells
    //     - overwrites what was previously an outer guard cell
    //
    //  Notes:
    //    - Simply copies the primitives of the last valid cell.
    //    - Assumes the parameter LogZoning flag is still valid
    // ============================================= //

    std::cout << "Extending grid" << std::endl;

    const unsigned int Nr_old = theDomain->Nr;
    if( Nr_old > 50000 )
    {
        throw std::runtime_error("too many zones");
    }

    theDomain->Nr *= 1.25;
    const unsigned int Nr_new = theDomain->Nr;

    theDomain->theCells = (struct cell *) realloc( theDomain->theCells , 
                                                   Nr_new*sizeof(struct cell) );

    int LogZoning;
    if ( this->trust_LogZoning_flag() == true)
    {
        LogZoning = theDomain->theParList.LogZoning;
    }
    else
    {
        const double tolerance = 1e-4;

        const double r_3 = theDomain->theCells[Nr_old-2].riph;
        const double r_2 = theDomain->theCells[Nr_old-3].riph;
        const double r_1 = theDomain->theCells[Nr_old-4].riph;


        const double dr_32 = r_3 - r_2;
        const double dr_21 = r_2 - r_1;

        if ( (std::abs(dr_32 - dr_21) / dr_21) < tolerance )
        {
            LogZoning = 0;
        }
        else if ( (std::abs((dr_32/r_3) - (dr_21/r_2)) / (dr_21/r_2))
                        < tolerance )
        {
            LogZoning = 1;
        }
        else
        {
            std::cerr << "Warning: extend_grid couldn't determine"
                      << " linear/log spacing; defaulting to linear"
                      << std::endl;
            LogZoning = 0;
        }
    }


    for( int i = Nr_old-2 ; i < Nr_new ; ++i  )
    {
        struct cell * c = &(theDomain->theCells[i]);
        const double rm =  theDomain->theCells[i-1].riph;
        const double dr_previous = rm - theDomain->theCells[i-2].riph;

        if( LogZoning == 0 )
        {
            c->riph = rm + dr_previous;
        }
        else 
        {
            c->riph = rm / (1 - (dr_previous / rm));
        }
    }

    calc_dr( theDomain );

    struct cell * c_old = &(theDomain->theCells[Nr_old-2]);
    for( int i = Nr_old-2 ; i < Nr_new ; ++i )
    {
        struct cell * c = &(theDomain->theCells[i]);
        c->wiph = c_old->wiph;

        for (int q=0 ; q < NUM_Q ; ++q)
        {
            c->prim[q] = c_old->prim[q];
        }

        const double rp = c->riph;
        const double rm = rp - c->dr;
        const double dV = get_dV( rp , rm );
        prim2cons( c->prim , c->cons , dV );
        cons2prim( c->cons , c->prim , dV );
        c->E_int_old = E_int_from_cons( c->cons );
        c->dV_old = dV;
    }    


    boundary( theDomain );

};

void Initial_Conditions::set_times( struct domain * theDomain )
{

    assert( std::is_sorted( theDomain->SNe.rbegin(),
                            theDomain->SNe.rend(),
                            sort_by_lifetime) );

    if ( theDomain->SNe.size() > 0 )
    {

        double t_first_SN = theDomain->SNe.back().lifetime;
        double t_last_SN  = theDomain->SNe.front().lifetime;

        theDomain->t      += t_first_SN;
        theDomain->t_init += t_first_SN;
        theDomain->t_fin  += t_last_SN;

    }
    else
    {
        // std::cerr << "Error: No SNe in this run. Exiting." << std::endl;
        // // no supernovae. For now, just kill the process
        // // but maybe figure out a better way to respond?
        // return 1; 
    }

}

void Initial_Conditions::set_output_prefix( struct domain * theDomain )
{

    uuid_t  id_binary;
    char id_ascii[36];
    uuid_generate(id_binary);
    uuid_unparse(id_binary, id_ascii);
    printf("generated uuid: %s \n", id_ascii);
    theDomain->output_prefix = std::string(id_ascii).append("_");

}


