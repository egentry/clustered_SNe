
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>

#include "structure.H"
#include "constants.H"
#include "misc.H"
#include "mass_loss.H"
#include "boundary.H" // boundary
#include "misc.H" // calc_prim

#include <boost/algorithm/string.hpp>    


Mass_Loss * select_mass_loss( std::string mass_loss_name )
{

    std::cout << "Using mass loss prescription: " << mass_loss_name << std::endl;

    boost::algorithm::to_lower(mass_loss_name);

    if ( mass_loss_name.compare(boost::to_lower_copy(No_Mass_Loss::class_name)) == 0 )
    {
        return new No_Mass_Loss;
    }

    if ( mass_loss_name.compare(boost::to_lower_copy(Uniform_Mass_Loss::class_name)) == 0 )
    {
        return new Uniform_Mass_Loss;
    }

    if ( mass_loss_name.compare(boost::to_lower_copy(Disappear_Mass_Loss::class_name)) == 0 )
    {
        return new Disappear_Mass_Loss;
    }

    throw std::invalid_argument("Mass loss name didn't match known names");

}

//********************** Mass_Loss (abstract base class) *********************//
Mass_Loss::Mass_Loss( const std::string name ) : name(name)
{}

double Mass_Loss::get_ejecta_mass( const double M_initial ) const
{
    // ============================================= //
    //
    //  Given an initial mass, find the ejecta mass (total)
    //  using the data of Heger & Woosley (2007)
    // 
    //  Inputs:
    //     - M_initial    - initial stellar mass [M_sol]
    //
    //  Outputs:
    //     - M_ejecta     - ejected mass [M_sol]
    //
    //  Side effects:
    //     - None
    //
    //  Notes:
    //     - This assumes that there *is* a supernova
    //     - (see Fig 8 of http://arxiv.org/pdf/1503.07522v1.pdf)
    //
    // ============================================= //

    const unsigned int grid_size = 32;
    const double M_initials[] = {  12.0,  13.0, 14.0,
                                   15.0,  16.0, 17.0,
                                   18.0,  19.0, 20.0,
                                   21.0,  22.0, 23.0,
                                   24.0,  25.0, 26.0,
                                   27.0,  28.0, 29.0,
                                   30.0,  31.0, 32.0,
                                   33.0,  35.0, 40.0,
                                   45.0,  50.0, 55.0,
                                   60.0,  70.0, 80.0,
                                  100.0, 120.0 };

    const double M_ejectas[] = {   9.389336,  10.244748, 10.973754,
                                  11.837161,  12.628000, 13.283933,
                                  13.973949,  14.867287, 15.229942,
                                  15.757695,  15.604063, 16.343219,
                                  16.144743,  15.389504, 14.995688,
                                  15.197527,  15.178074, 14.443274,
                                  14.193992,  14.117690, 13.952423,
                                  13.509647,  13.313397, 13.354588,
                                  15.109093,  10.222817,  8.271961,
                                   5.883033,   4.969768,  4.579023,
                                   4.227612,   4.386078 };
    if ( M_initial < M_initials[0] )
    {
        return M_initial - 1.5; // extrapolate Fig 3 of Woosley + Heger (2007)
    }

    for (unsigned int i=0 ; i<grid_size ; ++i)
    {
        if( M_initial < M_initials[i] )
        {
            // linearly interpolate
            return M_ejectas[i-1] 
                + (M_ejectas[i] - M_ejectas[i-1])  
                    * (M_initial     - M_initials[i-1]) 
                    / (M_initials[i] - M_initials[i-1]);
        }
    }

    return M_ejectas[grid_size-1];

}

double Mass_Loss::get_ejecta_mass_Z( const double M_initial ) const 
{
    // ============================================= //
    //
    //  Given an initial mass, find the ejecta mass (metals)
    //  using the data of Heger & Woosley (2007)
    // 
    //  Inputs:
    //     - M_initial    - initial stellar mass [M_sol]
    //
    //  Outputs:
    //     - M_ejecta_Z   - ejected mass of metals [M_sol]
    //
    //  Side effects:
    //     - None
    //
    //  Notes:
    //     - This assumes that there *is* a supernova
    //     - (see Fig 8 of http://arxiv.org/pdf/1503.07522v1.pdf)
    //
    // ============================================= //



    const unsigned int grid_size = 32;
    const double M_initials[] = {  12.0,  13.0, 14.0,
                                   15.0,  16.0, 17.0,
                                   18.0,  19.0, 20.0,
                                   21.0,  22.0, 23.0,
                                   24.0,  25.0, 26.0,
                                   27.0,  28.0, 29.0,
                                   30.0,  31.0, 32.0,
                                   33.0,  35.0, 40.0,
                                   45.0,  50.0, 55.0,
                                   60.0,  70.0, 80.0,
                                  100.0, 120.0 };

    const double M_ejecta_Zs[] = {  0.712336,  0.520048, 0.740254,
                                    1.026461,  1.282800, 1.549433,
                                    1.724449,  2.314587, 2.590042,
                                    2.916695,  3.095063, 3.709219,
                                    3.934743,  3.864504, 4.228688,
                                    4.699527,  5.152074, 5.776274,
                                    6.193992,  6.813690, 7.216423,
                                    7.386647,  8.174397, 9.856588,
                                   11.739093, 10.172817, 8.211961,
                                    5.813033,  4.829768, 4.409023,
                                    4.037612,  4.216078};


    if ( M_initial < M_initials[0] )
    {
        return M_ejecta_Zs[0];
    }

    for (unsigned int i=0 ; i<grid_size ; ++i)
    {
        if( M_initial < M_initials[i] )
        {
            // linearly interpolate
            return M_ejecta_Zs[i-1] 
                + (M_ejecta_Zs[i] - M_ejecta_Zs[i-1])  
                    * (M_initial     - M_initials[i-1]) 
                    / (M_initials[i] - M_initials[i-1]);
        } 
    }

    return M_ejecta_Zs[grid_size-1];

}


double Mass_Loss::get_wind_mass( const double M_initial ) const 
{
    // ============================================= //
    //
    //  Given an initial mass, find the total wind mass
    //  using the data of Heger & Woosley (2007)
    // 
    //  Inputs:
    //     - M_initial    - initial stellar mass [M_sol]
    //
    //  Outputs:
    //     - M_wind       - mass lost pre-supernova [M_sol]
    //
    //  Side effects:
    //     - None
    //
    //  Notes:
    //     - This assumes that there *is* a supernova
    //     - (see Fig 8 of http://arxiv.org/pdf/1503.07522v1.pdf)
    //
    // ============================================= //



    const unsigned int grid_size = 32;
    const double M_initials[] = {  12.0,  13.0, 14.0,
                                   15.0,  16.0, 17.0,
                                   18.0,  19.0, 20.0,
                                   21.0,  22.0, 23.0,
                                   24.0,  25.0, 26.0,
                                   27.0,  28.0, 29.0,
                                   30.0,  31.0, 32.0,
                                   33.0,  35.0, 40.0,
                                   45.0,  50.0, 55.0,
                                   60.0,  70.0, 80.0,
                                  100.0, 120.0 };

    const double M_winds[] = {   1.086375,   1.238961,  1.486873,
                                 1.632861,   1.816102,  2.085835,
                                 2.258426,   2.599328,  3.188975,
                                 3.616268,   4.586832,  5.144166,
                                 6.223336,   7.510593,  8.849613,
                                 9.776460,  10.875817, 12.867450,
                                14.131232,  15.385821, 16.551013,
                                17.822844,  19.897414, 24.516315,
                                27.569192,  37.675491, 45.046299,
                                52.009352,  63.382721, 73.839761,
                                93.962149, 114.005312 };


    if ( M_initial < M_initials[0] )
    {
        // it'll be less than 1 solar mass, so whatever we put
        // probably won't matter
        return 0;
    }

    for (unsigned int i=0 ; i<grid_size ; ++i)
    {
        if( M_initial < M_initials[i] )
        {
            // linearly interpolate
            return M_winds[i-1] 
                + (M_winds[i] - M_winds[i-1])  
                    * (M_initial     - M_initials[i-1]) 
                    / (M_initials[i] - M_initials[i-1]);
        } 
    }

    return M_winds[grid_size-1];

}

//********************** Uniform_Mass_Loss *********************//
const std::string Uniform_Mass_Loss::class_name = "uniform";

Uniform_Mass_Loss::Uniform_Mass_Loss() : Mass_Loss( class_name )
{}

void Uniform_Mass_Loss::add_mass_loss( struct domain * theDomain , 
                                       const double dt ) const 
{
    // ============================================= //
    //
    //  Simply assume a constant mass-loss rate,
    //   at a constant wind velocity, 
    //   with the metallicity being the raw metallicity of the star
    // 
    //  Inputs: 
    //     - theDomain      - the standard domain struct used throughout
    //     - dt             - dt of the current timestep
    //
    //  Outputs:
    //     - None
    //
    //  Side effects:
    //     - Adds mass, metals, momentum, energy to the cons
    //
    //  Notes:
    //      - Overall mass loss should have a relatively small effect,
    //        but for the most massive SNe, pre-SNe mass loss 
    //        (e.g. during Wolf Rayet phase) can change the SNR dynamics
    //      - Only adds mechanical energy
    //          - I should probably add thermal energy?
    //
    // ============================================= //

    if( theDomain->rank != 0) return;

    const double wind_velocity = 1000 * (1000 * 100); // 1000 km / s


    const int n_guard_cell = 1;
    struct cell * blast_cell = &(theDomain->theCells[n_guard_cell]);


    for ( unsigned int i=0 ; i <theDomain->SNe.size() ; ++i )
    {
        supernova SN = theDomain->SNe[i];
        const double mass_loss = SN.mass_winds *  (dt / SN.lifetime); 

        blast_cell->cons[DDD] += mass_loss;
        blast_cell->cons[SRR] += mass_loss * wind_velocity;

        const double T_wind = 1e4; // K
        const double mu = get_mean_molecular_weight( theDomain->metallicity );
        blast_cell->cons[TAU] += mass_loss * .5 * std::pow(wind_velocity, 2) 
                              +  mass_loss * k_boltzmann * T_wind / (mu * m_proton);
        blast_cell->cons[ZZZ] += mass_loss * theDomain->metallicity;

    }

    // now we need to background within substep()
    // since we changed the cons, but not the prims
    calc_prim( theDomain );
    boundary(  theDomain );

    return;
}

//********************** No_Mass_Loss *********************//
const std::string No_Mass_Loss::class_name = "none";


No_Mass_Loss::No_Mass_Loss() : Mass_Loss( class_name )
{}

double No_Mass_Loss::get_wind_mass( const double M_initial ) const
{
    return 0.0;
}

double No_Mass_Loss::get_ejecta_mass( const double M_initial ) const
{
    // approximated from Fig 3 of Woosley + Heger (2007)
    // this is also consistent with the Agora simulations
    const double M_remnant = 1.4; // solar masses, 


    return M_initial - M_remnant;
}

double No_Mass_Loss::get_ejecta_mass_Z( const double M_initial ) const
{
    // See Section 3.5 of the 2013 Agora paper: arXiv:1308.2669

    const double M_Fe = 0.375 * std::exp(-17.94 / M_initial); // solar masses
    const double M_O  = 27.66 * std::exp(-51.81 / M_initial); // solar masses

    const double M_ejecta_Z = (2.09 * M_O) + (1.06 * M_Fe);

    return M_ejecta_Z;

}

//********************** Disappear_Mass_Loss *********************//
const std::string Disappear_Mass_Loss::class_name = "disappear";

Disappear_Mass_Loss::Disappear_Mass_Loss() : Mass_Loss( class_name )
{}

double Disappear_Mass_Loss::get_wind_mass( const double M_initial ) const
{
    return 0.0;
}



