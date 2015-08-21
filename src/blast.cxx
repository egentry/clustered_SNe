
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
#include "constants.H"

#include "blast.H"

using namespace boost;
using namespace boost::filesystem;

double get_dV( double , double );
void calc_prim( struct domain * );
void boundary( struct domain * );
int add_single_blast( struct domain * theDomain , const double E_blast )
{

    // ============================================= //
    //
    //  Adds a blast wave of a given energy
    //  to the innermost valid cell
    // 
    //  Inputs:
    //     - theDomain    - the standard domain struct used throughout
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

    const double M_blast = 3 * M_sun;

    // code is not yet set for spreading SNe over multiple cells
    // would need to determine how to split the energy correctly

    struct cell * c = &(theDomain->theCells[n_guard_cell]);

    double previous_metallicity = c->cons[ZZZ] / c->cons[DDD];

    c->cons[DDD] += M_blast;
    c->cons[TAU] += E_blast;

    // for now assume that the metallicity of the ejecta is the same
    // as the background metallicity
    // later we'll want to actually inject metals
    c->cons[ZZZ] += M_blast * previous_metallicity;


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

    assert( std::is_sorted( theDomain->SNe_times.rbegin(),
                            theDomain->SNe_times.rend()   ) );

    unsigned int num_SNe = 0;

    while ( (theDomain->SNe_times.size() > 0) && 
            (theDomain->SNe_times.back() < theDomain->t) )
    {
        // double SNe_time = theDomain->SNe_times.back();
        theDomain->SNe_times.pop_back();

        ++num_SNe;
    }

    int error = 0;
    for (int i = 0; i<num_SNe; ++i)
    {
        error += add_single_blast( theDomain );
    }
    if (num_SNe > 0)
        std::cout << num_SNe << " SN(e) just went off" << std::endl;

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


int get_SNe( const double cluster_mass, std::vector<double>&SNe_times,
                      const double metallicity,
                      const unsigned int seed)
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
    //     - SNe_times      - vector to contain 
    //     - metallicity    - metallicity mass fraction (e.g. solar = 0.02)
    //                      - Default = 0.02
    //     - seed           - seeds the stochasticity;
    //
    //  Outputs:
    //     - error          - 0 for successful execution
    //                      - 1 for failure (should cause this run to quietly stop)
    //
    //  Side effects:
    //     - overwrites SNe_times completely
    //
    //  Notes:
    //
    // ============================================= //

    SNe_times.resize(0);

    typedef boost::random::mt19937 rng_type; 
    rng_type *rng;

    rng = new rng_type(seed); 

    // SLUG_DIR should be macro'd in at compile time
    boost::filesystem::path slug_path(SLUG_DIR);

    boost::filesystem::path imf_path("lib/imf/kroupa.imf");
    imf_path = slug_path / imf_path;

    boost::filesystem::path tracks_path("lib/tracks/Z0140v00.txt");
    tracks_path = slug_path / tracks_path;

    // while( SNe_times.size() != 2)
    // {
        slug_PDF imf(imf_path.string().c_str(), rng, true); 
        std::vector<double> star_masses; 
        imf.drawPopulation(cluster_mass, star_masses);

        slug_tracks tracks(tracks_path.string().c_str(), metallicity); 
        for (auto mass: star_masses)
        {
            if ( mass > 8.0 ) // must be 8 solar masses to be consistent with slug
            {
                SNe_times.push_back(tracks.star_lifetime(mass)); // [yr]
            }
        }
    // }
    // sort in reverse order, so we can pop off SNe as they happen
    std::sort(SNe_times.rbegin(), SNe_times.rend());

    for (double& SNe_time : SNe_times )
    {
        SNe_time *= yr; // convert to [s]
    }

    std::cout << "Num SNe: " << SNe_times.size() << std::endl;
    std::cout << "Cluster mass: " << cluster_mass << std::endl;

    std::cout << "SNe times: " << std::endl;
    for (auto t : SNe_times)
        std::cout << t << "[s]" << std::endl;

    return 0;
}