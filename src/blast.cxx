
#include <iostream>
#include <string>
#include <random>
#include <assert.h>


#include "slug_PDF.H"
#include "slug_tracks.H"
#include <vector>
#include <algorithm>    // std::sort
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>


#include <assert.h>

#include "constants.H" // pc, m_proton, k_boltzmann, yr
#include "structure.H"
#include "boundary.H"
#include "constants.H"
#include "geometry.H" // get_dV
#include "misc.H" // calc_prim
#include "mass_loss.H" // Mass_Loss, get_ejecta_mass, etc

#include "blast.H"

using namespace boost;
using namespace boost::filesystem;

double calc_R_E(const double rho, const double Z, 
           const double P_0, const double E_blast)
{
    // ============================================= //
    //
    //  Determines R_E of Stinson et al. 2006, equation 9
    // 
    //  Inputs:
    //     - rho          - Ambient mas density within R_E [g cm**-3]
    //     - Z            - Ambient gas-phase metallicity
    //     - P_0          - Ambient [volume-averaged] pressure within R_E [erg cm**-3]
    //     - E_blast      - the energy to be injected [ergs]
    //
    //  Outputs:
    //     - R_E          - radius of maximum extent [cm]
    //
    //  Side effects:
    //     None
    //
    // ============================================= //

    const double Y = .23; // helium fraction
    const double X = 1 - Y - Z; // hydrogen mass fraction

    const double n_0 = rho * X / m_proton; // * HYDROGEN number density [cm**-3]
    const double P_04 = 1e-4 * P_0 / k_boltzmann; 
    const double E_51 = E_blast / 1e51;

    // Stinson et al 2006, Equation 9 (from Chevalier 1974)
    const double R_E =  std::pow(10, 1.74) 
        * std::pow(E_51,  0.32)
        * std::pow(n_0,  -0.16)
        * std::pow(P_04, -0.20)
        * pc;

    return R_E;

}

double calc_t_E(const double rho, const double Z, 
           const double P_0, const double E_blast)
{
    // ============================================= //
    //
    //  Determines t_E of Stinson et al. 2006, equation 10
    // 
    //  Inputs:
    //     - rho          - Ambient mas density within R_E [g cm**-3]
    //     - Z            - Ambient gas-phase metallicity
    //     - P_0          - Ambient [volume-averaged] pressure within R_E [erg cm**-3]
    //     - E_blast      - the energy to be injected [ergs]
    //
    //  Outputs:
    //     - t_E          - time of maximum extent [s]
    //
    //  Side effects:
    //     None
    //
    //  Notes:
    //     - Corresponds to the time when the SN should get to R_E
    // ============================================= //

    const double Y = .23; // helium fraction
    const double X = 1 - Y - Z; // hydrogen mass fraction

    const double n_0 = rho * X / m_proton; // * HYDROGEN number density [cm**-3]
    const double P_04 = 1e-4 * P_0 / k_boltzmann; 
    const double E_51 = E_blast / 1e51;

    // Stinson et al 2006, Equation 10 (from Chevalier 1974)
    const double t_E =  std::pow(10, 5.92) 
        * std::pow(E_51,  0.31)
        * std::pow(n_0,   0.27)
        * std::pow(P_04, -0.64)
        * yr;

    // printf("E_51 = %f \n", E_51);
    // printf("n_0 = %f \n", n_0);
    // printf("P_04 = %f \n", P_04);
    // printf("t_E = %e \n \n", t_E / yr);

    return t_E;

}

double calc_t_max(const double rho, const double Z, 
                  const double P_0, const double E_blast)
{
    // ============================================= //
    //
    //  Determines t_max of Stinson et al. 2006, equation 11
    // 
    //  Inputs:
    //     - rho          - Ambient mas density within R_E [g cm**-3]
    //     - Z            - Ambient gas-phase metallicity
    //     - P_0          - Ambient [volume-averaged] pressure within R_E [erg cm**-3]
    //     - E_blast      - the energy to be injected [ergs]
    //
    //  Outputs:
    //     - t_max          - lifetime of low-density shell [s]
    //
    //  Side effects:
    //     None
    //
    //  Notes:
    //     - Corresponds to the time when the SN should get to R_E
    // ============================================= //

    const double Y = .23; // helium fraction
    const double X = 1 - Y - Z; // hydrogen mass fraction

    const double n_0 = rho * X / m_proton; // * HYDROGEN number density [cm**-3]
    const double P_04 = 1e-4 * P_0 / k_boltzmann; 
    const double E_51 = E_blast / 1e51;

    // Stinson et al 2006, Equation 10 (from Chevalier 1974)
    const double t_max =  std::pow(10, 6.85) 
        * std::pow(E_51,  0.32)
        * std::pow(n_0,   0.34)
        * std::pow(P_04, -0.70)
        * yr;

    return t_max;

}

double guassian(const double r, const double sigma)
{
    return std::exp(- .5 * std::pow(r / sigma, 2.));
}



int add_single_blast( struct domain * theDomain , const double M_blast ,
                                                  const double M_blast_Z,
                                                  const double E_blast )
{

    // ============================================= //
    //
    //  Adds a blast wave of a given energy
    //  to the innermost valid cell
    // 
    //  Inputs:
    //     - theDomain    - the standard domain struct used throughout
    //     - M_blast      - total mass ejected [g]
    //     - M_blast_Z    - mass of metals ejected [g] 
    //     - E_blast      - the energy to be injected [ergs]
    //                    - default = 1e51 ergs
    //
    //  Outputs:
    //     - error        - 0 for successful execution
    //                    - 1 for failure (should cause this rank to quietly stop)
    //
    //  Side effects:
    //     - overwrites the cons AND prims
    //       of the innermost non-guard cell
    //     - enforces boundary condition
    //
    //  Notes:
    //     - Should only be used before full timesteps
    //
    // ============================================= //

    const int n_guard_cell = 1;

    // iteratively solve to find the correct R_E
    double Mass = 0;
    double Mass_Z = 0;
    double Volume = 0;
    double Pressure_times_Volume = 0;

    double t_E = 0;
    double R_E = 0;
    for( int i=theDomain->Ng ; i<theDomain->Nr-theDomain->Ng ; ++i )
    {
        struct cell * c = &(theDomain->theCells[i]);

        Mass += c->cons[DDD];
        Mass_Z += c->cons[ZZZ];

        const double dV = get_dV(c->riph, c->riph - c->dr);

        Volume += dV;

        // check: do I need to re-compute the primitives?
        Pressure_times_Volume += c->prim[PPP] * dV;


        R_E = calc_R_E(Mass / Volume, Mass_Z / Mass,
                       Pressure_times_Volume / Volume, E_blast);
        t_E = calc_t_E(Mass / Volume, Mass_Z / Mass,
                       Pressure_times_Volume / Volume, E_blast);

        if (R_E < c->riph)
        {
            printf("breaking at i= %d \n",i);

            break;
        }

    }

    printf("R_E solution found!\n");
    printf("R_E    = %e pc \n", R_E / pc);
    printf("t_E    = %e yr \n", t_E / yr);

    printf("Mass   = %e M_solar \n", Mass / M_sun);
    printf("Mass_Z = %e M_solar \n", Mass_Z / M_sun);
    printf("rho   = %e g cm**-3 \n", Mass / Volume);
    printf("Z      = %e \n", Mass_Z / Mass);
    printf("Volume = %e pc**3 \n", Volume / std::pow(pc, 3.));
    printf("Pressure / k_B = %e K\n", Pressure_times_Volume / Volume / k_boltzmann);


    assert(R_E > 0);

    // calculate the array of weights


    const double Mass_kernel = 3e5 * M_sun; // ideally change this to be set at run-time using first SN
    double Mass_tmp = 0;
    double R_kernel = 0;
    for( int i=theDomain->Ng ; i<theDomain->Nr-theDomain->Ng ; ++i )
    {
        struct cell * c = &(theDomain->theCells[i]);

        R_kernel = c->riph;
        Mass_tmp += c->cons[DDD];

        if(Mass_tmp > Mass_kernel)
        {
            printf("Injection cutoff: %f pc \n", R_kernel / pc);
            break;
        }
    }


    const double sigma = .1 * R_kernel;

    double weights[theDomain->Nr];
    for( int i=0; i<theDomain->Nr; ++i ) weights[i] = 0;


    double total_weight = 0;
    for( int i=theDomain->Ng ; i<theDomain->Nr-theDomain->Ng ; ++i )
    {
        struct cell * c = &(theDomain->theCells[i]);

        const double r = c->riph;

        if (r <= R_kernel)
        {
            weights[i] = guassian(r, sigma) * c->cons[DDD];
        }
        total_weight += weights[i];
    }

    // normalize the weights

    for( int i=theDomain->Ng ; i<theDomain->Nr-theDomain->Ng ; ++i )
    {
        weights[i] /= total_weight;
    }

    // add the energy, mass and metals using the weights
    for( int i=theDomain->Ng ; i<theDomain->Nr-theDomain->Ng ; ++i )
    {
        struct cell * c = &(theDomain->theCells[i]);

        c->cons[DDD] += weights[i] * M_blast;
        c->cons[TAU] += weights[i] * E_blast;
        c->cons[ZZZ] += weights[i] * M_blast_Z;
    }

    // shutoff cooling for applicable cells
    const double R_cutoff = std::min(R_E, R_kernel);
    printf("R_E      = %f pc \n", R_E / pc);
    printf("R_kernel = %f pc \n", R_kernel / pc);
    printf("R_cutoff = %f pc \n", R_cutoff / pc);

    for( int i=theDomain->Ng ; i<theDomain->Nr-theDomain->Ng ; ++i )
    {
        struct cell * c = &(theDomain->theCells[i]);

        if (c->riph < R_cutoff)
        {
            c->cooling_active = 0;
            c->shutoff_cooling_until = std::max(c->shutoff_cooling_until, 
                                                t_E + theDomain->t);
        }        

    }



    // code is not yet set for spreading SNe over multiple cells
    // would need to determine how to split the energy correctly

    // struct cell * c = &(theDomain->theCells[n_guard_cell]);

    // c->cons[DDD] += M_blast;
    // c->cons[TAU] += E_blast;

    // // for now assume that the metallicity of the ejecta is the same
    // // as the background metallicity
    // // later we'll want to actually inject metals
    // c->cons[ZZZ] += M_blast_Z;


    // now we need to background within substep()
    // since we changed the cons, but not the prims
    calc_prim( theDomain );
    boundary( theDomain );


    return 0;
}

int add_blasts( struct domain * theDomain )
{
    // ============================================= //
    //
    //  Adds as many blast waves as necessary within a timestep
    //    (doesn't guarantee that these blasts should have happened
    //     within a given dt, just that all the blasts that should have
    //     gone off up to the current time have now gone off)
    // 
    //  Inputs:
    //     - theDomain    - the standard domain struct used throughout
    //
    //  Outputs:
    //     - error        - 0 for successful execution
    //                    - 1 for failure (should cause this rank to quietly stop)
    //
    //  Side effects:
    //    If there is 1 or more supernovae:
    //     - pop
    //     - also incurs side effects of add_single_blast
    //
    //  Notes:
    //     - Should only be used before full timesteps,
    //       else you might run into CFL errors
    //
    // ============================================= //

    const double E_SN = 1e51; // less than 1e51 because we turnoff cooling

    assert( std::is_sorted( theDomain->SNe.rbegin(),
                            theDomain->SNe.rend(),
                            sort_by_lifetime) );

    unsigned int num_SNe = 0;
    int error = 0;
    while ( (theDomain->SNe.size() > 0) && 
            (theDomain->SNe.back().lifetime < theDomain->t) )
    {
        supernova tmp = theDomain->SNe.back();
        error += add_single_blast( theDomain, 
                                   tmp.mass_ejecta, tmp.mass_ejecta_Z,
                                   E_SN);

        theDomain->SNe.pop_back();
        ++num_SNe;
    }

    if ( num_SNe > 0 )
    {
        std::cout << num_SNe << " SN(e) just went off" << std::endl << std::flush;
    }

    if ( error > 0 )
    {
        std::cerr << "Error in adding blasts to the simulation" << std::endl;
        return 1;
    }
    else
    {
        return 0;
    }
}


std::vector<supernova> get_SNe( const double cluster_mass ,
                                const double metallicity ,
                                const unsigned int seed ,
                                const Mass_Loss * mass_loss )
{
    // ============================================= //
    //
    //  Gets all the necessary information about the stochastic supernovae.
    //    This information includes:
    //      - List of times for each supernova after cluster formation
    //    This information does not include:
    //      - Energy of each supernova
    //
    //  Inputs:
    //     - cluster_mass   - cluster mass in solar masses
    //     - metallicity    - metallicity mass fraction (e.g. solar = 0.02)
    //                      - Default = 0.02
    //     - seed           - seeds the stochasticity;
    //     - mass_loss      - sets the wind mass and ejecta mass
    //
    //  Outputs:
    //     - SNe            - contains relevant information of the blasts
    //
    //  Side effects:
    //     - None
    //
    //  Notes:
    //     - Assumes a Kroupa IMF
    //     - Uses metallicity dependent stellar evolution tracks
    //     - Assumes a SNe if and only if a star is over 8 solar masses
    //          - not all of these stars will explode in actuality
    //     - Assumes solar metallicity progenitors for ejecta mass
    //
    // ============================================= //

    typedef boost::random::mt19937 rng_type; 
    rng_type *rng;

    rng = new rng_type(seed); 

    // SLUG_DIR should be macro'd in at compile time
    boost::filesystem::path slug_path(SLUG_DIR);

    boost::filesystem::path imf_path("lib/imf/kroupa.imf");
    imf_path = slug_path / imf_path;

    boost::filesystem::path tracks_path("lib/tracks/Z0140v00.txt");
    tracks_path = slug_path / tracks_path;

    std::vector<double> star_masses; 
    slug_PDF imf(imf_path.string().c_str(), rng, true); 
    imf.drawPopulation(cluster_mass, star_masses);

    slug_tracks tracks(tracks_path.string().c_str(), metallicity); 
    std::vector<supernova> SNe;
    for (auto mass: star_masses)
    {
        if ( mass > 8.0 ) // must be 8 solar masses to be consistent with slug
        {
            supernova tmp;
            tmp.mass          = mass * M_sun;                    // [g]
            tmp.lifetime      = tracks.star_lifetime(mass) * yr; // [s]
            tmp.mass_ejecta   = mass_loss->get_ejecta_mass(  mass) * M_sun; // [g]
            tmp.mass_ejecta_Z = mass_loss->get_ejecta_mass_Z(mass) * M_sun; // [g]
            tmp.mass_winds    = mass_loss->get_wind_mass(mass) * M_sun; // [g]
            SNe.push_back(tmp);
        }
    }

    // sort in reverse order, so we can pop off SNe as they happen
    std::sort(SNe.rbegin(), SNe.rend(), sort_by_lifetime);

    std::cout << "Num SNe: " << SNe.size() << std::endl;
    std::cout << "Cluster mass: " << cluster_mass  << " M_sun" << std::endl;

    return SNe;
}


bool sort_by_lifetime( const supernova &a, const supernova &b)
{
    return (a.lifetime < b.lifetime);
}

