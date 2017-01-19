
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

#include "Hydro/euler.H" // E_int_from_cons

#include "blast.H"

using namespace boost;
using namespace boost::filesystem;


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

    // code is not yet set for spreading SNe over multiple cells
    // would need to determine how to split the energy correctly

    struct cell * c = &(theDomain->theCells[n_guard_cell]);

    if (c->multiphase==0)
    {
        c->cons_cold[DDD] = c->cons[DDD];
        c->cons_cold[TAU] = c->cons[TAU];
        c->cons_cold[ZZZ] = c->cons[ZZZ];
    }

    c->multiphase = 1;
    c->cons_hot[DDD] += M_blast;
    c->cons_hot[TAU] += E_blast;
    c->cons_hot[ZZZ] += M_blast_Z;

    c->cons[DDD] += M_blast;
    c->cons[TAU] += E_blast;
    // for now assume that the metallicity of the ejecta is the same
    // as the background metallicity
    // later we'll want to actually inject metals
    c->cons[ZZZ] += M_blast_Z;

    c->x_hot  = c->cons_hot[DDD]  / c->cons[DDD];
    c->x_cold = c->cons_cold[DDD] / c->cons[DDD];

    c->y_hot  = E_int_from_cons(c->cons_hot)  / E_int_from_cons(c->cons);
    c->y_cold = E_int_from_cons(c->cons_cold) / E_int_from_cons(c->cons);

    c->z_hot  = c->cons_hot[ZZZ]  / c->cons_hot[DDD];
    c->z_cold = c->cons_cold[ZZZ] / c->cons_cold[DDD];

    c->E_int_initial = E_int_from_cons(c->cons);
    c->E_kin_initial = c->cons[TAU] - c->E_int_initial;


    if (c->multiphase)
    {
        const double rel_tol = 1e-5; // relative tolerance for float comparisons
        if (c->cons_hot[DDD] < 0)
        {
            printf("------ERROR in add_single_blast()------- \n");
            printf("hot gas mass less than 0\n");
            printf("c->cons_hot[DDD] = %e \n", c->cons_hot[DDD]);
            assert( c->cons_hot[DDD] > 0 );
        }

        if (c->cons_cold[DDD] < 0)
        {
            printf("------ERROR in add_single_blast()------- \n");
            printf("cold gas mass less than 0\n");
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            assert( c->cons_cold[DDD] > 0 );
        }

        if ( std::abs(1-( (c->cons_cold[DDD] + c->cons_hot[DDD])/c->cons[DDD])) > rel_tol)
        {
            printf("------ERROR in add_single_blast()------- \n");
            printf("cold mass + hot mass =/= total mass\n");
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);
            assert(  std::abs(1-( (c->cons_cold[DDD] + c->cons_hot[DDD])/c->cons[DDD])) <= rel_tol);
        }

        if (c->cons_hot[TAU] < 0)
        {
            printf("------ERROR in add_single_blast()------- \n");
            printf("hot gas energy less than 0\n");
            printf("c->cons_hot[TAU] = %e \n", c->cons_hot[TAU]);
            assert( c->cons_hot[TAU] > 0 );
        }

        if (c->cons_cold[TAU] < 0)
        {
            printf("------ERROR in add_single_blast()------- \n");
            printf("cold gas energy less than 0\n");
            printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
            assert( c->cons_cold[TAU] > 0 );
        }

        if ( std::abs(1-( (c->cons_cold[TAU] + c->cons_hot[TAU])/c->cons[TAU])) > rel_tol)
        {
            printf("------ERROR in add_single_blast()------- \n");
            printf("cold energy + hot energy =/= total energy\n");
            printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
            printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
            printf("c->cons[TAU]      = %e \n", c->cons[TAU]);
            assert(  std::abs(1-( (c->cons_cold[TAU] + c->cons_hot[TAU])/c->cons[TAU])) <= rel_tol);
        }
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
    int error = 0;
    while ( (theDomain->SNe.size() > 0) && 
            (theDomain->SNe.back().lifetime < theDomain->t) )
    {
        supernova tmp = theDomain->SNe.back();
        error += add_single_blast( theDomain, 
                                   tmp.mass_ejecta, tmp.mass_ejecta_Z);

        theDomain->SNe.pop_back();
        ++num_SNe;
    }

    if ( num_SNe > 0 )
    {
        std::cout << num_SNe << " SN(e) just went off" << std::endl;
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

