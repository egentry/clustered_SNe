
#include <iostream>
#include <string>
#include <random>

#include "slug_PDF.H"
#include "slug_tracks.H"
#include <vector>
#include <algorithm>    // std::sort
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>


#include <assert.h>

#include "structure.H"
#include "boundary.H"
#include "constants.H"
#include "geometry.H" // calc_dV
#include "misc.H" // calc_prim
#include "mass_loss.H" // Mass_Loss, get_ejecta_mass, etc

#include "blast.H"

using namespace boost;
using namespace boost::filesystem;


int add_single_blast( struct domain * theDomain , const double M_blast ,
                                                  const double M_blast_Z ,
                                                  const double E_blast ,
                                                  const int N_ngb )
{

    // ============================================= //
    //
    //  Adds a blast wave of a given mass, metals and energy
    //  to the innermost N_ngb cells (split evenly regardless
    //  of existing mass or volume).
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
    //       of the innermost N_ngb non-guard cells
    //     - enforces boundary condition
    //
    //  Notes:
    //     - Should only be used before full timesteps
    //
    // ============================================= //

    const int n_guard_cell = theDomain->Ng;

    // code is not yet set for spreading SNe over multiple cells
    // would need to determine how to split the energy correctly

    for(int i = n_guard_cell ; i < n_guard_cell + N_ngb ; ++i )
    {
        struct cell * c = &(theDomain->theCells[i]);

        c->cons[DDD] += M_blast / N_ngb;
        c->cons[TAU] += E_blast / N_ngb;

        // for now assume that the metallicity of the ejecta is the same
        // as the background metallicity
        // later we'll want to actually inject metals
        c->cons[ZZZ] += M_blast_Z / N_ngb;
    }



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

    assert( std::is_sorted( theDomain->SNe.rbegin(),
                            theDomain->SNe.rend(),
                            sort_by_lifetime) );

    unsigned int num_SNe = 0;
    const double E_SN = 1e51;
    const double N_ngb = 3;
    assert(N_ngb==3); // I'm hard coding the RNG results via input file
    double E_blasts = 0;
    int error = 0;
    while ( (theDomain->SNe.size() > 0) && 
            (theDomain->SNe.back().lifetime < theDomain->t) )
    {
        supernova tmp = theDomain->SNe.back();
        E_blasts += E_SN; // add later
        // just add mass and metals for now
        error += add_single_blast( theDomain, 
                                   tmp.mass_ejecta, tmp.mass_ejecta_Z,
                                   0 , N_ngb );

        theDomain->SNe.pop_back();
        ++num_SNe;
    }

    if( num_SNe == 0 )
    {
        return 0;
    }

    // now deal with stochastic injection
    const double dT_blast = theDomain->theParList.dT_blast;
    const double gamma = theDomain->theParList.Adiabatic_Index;
    const double mu = 2; // assume fully ionized, purely metals
    const double de = k_boltzmann * dT_blast / ((gamma-1) * mu * m_proton);

    double m_kernel = 0;
    for(int i = theDomain->Ng; i < theDomain->Ng + N_ngb ; ++i)
    {
        const struct cell * c = &(theDomain->theCells[i]);
        std::cout << "mass[i=" << i << "] = " << c->cons[DDD] << std::endl;
        m_kernel += c->cons[DDD];
    }

    const double p = E_blasts / de / m_kernel;
    std::cout << "p = " << p << std::endl;

    double dE_expected = 0;
    for(int i = theDomain->Ng; i < theDomain->Ng + N_ngb ; ++i)
    {
        const struct cell * c = &(theDomain->theCells[i]);
        std::cout << "p*dE[i=" << i << "] = " << de*p*c->cons[DDD] << std::endl;
        dE_expected += de*p*c->cons[DDD];
    }
    std::cout << "dE_expected = " << dE_expected << std::endl;

    // now stochastically add energy
    double dE_actually_added = 0;
    // std::mt19937 gen(theDomain->theParList.stochastic_seed);
    // std::uniform_real_distribution<> unif(0.0, 1.0);

    std::cout << "Add energy to respective cells: ";
    std::cout << theDomain->theParList.add_energy_cell_1 << " ";
    std::cout << theDomain->theParList.add_energy_cell_2 << " ";
    std::cout << theDomain->theParList.add_energy_cell_3 << std::endl;

    std::cout << "Likelihood: " << std::pow(p, 
        theDomain->theParList.add_energy_cell_1  
        + theDomain->theParList.add_energy_cell_2 
        + theDomain->theParList.add_energy_cell_3) 
    * std::pow(1-p, 
        N_ngb - (theDomain->theParList.add_energy_cell_1  
        + theDomain->theParList.add_energy_cell_2 
        + theDomain->theParList.add_energy_cell_3))
         << std::endl;

    if(theDomain->theParList.add_energy_cell_1)
    {
        const int i = theDomain->Ng + 0;
        struct cell * c = &(theDomain->theCells[i]);

        std::cout << "adding energy to cell i=" << i << std::endl;
        std::cout << "Energy before: " << c->cons[TAU] << std::endl;
        c->cons[TAU] += de * c->cons[DDD];
        dE_actually_added += de * c->cons[DDD];
        std::cout << "Energy after : " << c->cons[TAU] << std::endl;
    }
    if(theDomain->theParList.add_energy_cell_2)
    {
        const int i = theDomain->Ng + 1;
        struct cell * c = &(theDomain->theCells[i]);

        std::cout << "adding energy to cell i=" << i << std::endl;
        std::cout << "Energy before: " << c->cons[TAU] << std::endl;
        c->cons[TAU] += de * c->cons[DDD];
        dE_actually_added += de * c->cons[DDD];
        std::cout << "Energy after : " << c->cons[TAU] << std::endl;
    }
    if(theDomain->theParList.add_energy_cell_3)
    {
        const int i = theDomain->Ng + 2;
        struct cell * c = &(theDomain->theCells[i]);

        std::cout << "adding energy to cell i=" << i << std::endl;
        std::cout << "Energy before: " << c->cons[TAU] << std::endl;
        c->cons[TAU] += de * c->cons[DDD];
        dE_actually_added += de * c->cons[DDD];
        std::cout << "Energy after : " << c->cons[TAU] << std::endl;
    }

    // for(int i = theDomain->Ng; i < theDomain->Ng + N_ngb ; ++i)
    // {
    //     // const double r = unif(gen);
    //     std::cout << "cell i=" << i << " r = " << r << " p = " << p << std::endl;
    //     struct cell * c = &(theDomain->theCells[i]);

    //     if(r <= p)
    //     {
    //         std::cout << "adding energy to cell i=" << i << std::endl;
    //         std::cout << "Energy before: " << c->cons[TAU] << std::endl;
    //         c->cons[TAU] += de * c->cons[DDD];
    //         dE_actually_added += de * c->cons[DDD];
    //         std::cout << "Energy after : " << c->cons[TAU] << std::endl;
    //     }
    // }


    // now we need to background within substep()
    // since we changed the cons, but not the prims
    calc_prim( theDomain );
    boundary( theDomain );



    double fractional_E_var = 0;
    for(int i = theDomain->Ng; i < theDomain->Ng + N_ngb ; ++i)
    {
        const struct cell * c = &(theDomain->theCells[i]);
        fractional_E_var += std::pow(c->cons[DDD], 2);
    }
    fractional_E_var *= ((m_kernel * de/E_blasts) - 1) / std::pow(m_kernel, 2);

    if ( num_SNe > 0 )
    {
        std::cout << num_SNe << " SN(e) just went off" << std::endl;
        std::cout << "Using dT_blast = " << dT_blast << " K" << std::endl;
        std::cout << "Using de = " << de << " erg g^-1" << std::endl;
        std::cout << "m_kernel = " << m_kernel / M_sun << " M_sun" << std::endl;
        std::cout << "m_kernel * de = " << m_kernel * de << " erg" << std::endl;

        std::cout << "dE goal   : " << E_blasts << " erg" << std::endl;
        std::cout << "dE actual : " << dE_actually_added << " erg" << std::endl;
        std::cout << "fractional energy variance: " << fractional_E_var << std::endl;
        std::cout << "probability that at least 1 cell gets energy: " << 1 - std::pow(1-p, N_ngb) << std::endl;
        std::cout << std::flush;
    }

    return 1;

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

