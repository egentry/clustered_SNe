
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


    const double p_blast = std::sqrt(2 * E_blast * M_blast);

    // Calculate weights of cells
    const int N_neighbors = std::round(std::pow(32, 1./3.));

    const double weight = 1. / N_neighbors;

    const double dM = weight * M_blast;
    const double dp_prime = weight * p_blast;
    const double dE = weight * E_blast;
    const double dMZ = weight * M_blast_Z;

    double cons_added[NUM_Q];
    for (int q=0; q<NUM_Q ; ++q) cons_added[q] = 0;

    for (int i=0 ; i<N_neighbors ; ++i)
    {
        struct cell * c = &(theDomain->theCells[n_guard_cell + i]);

        const double E_int_old = E_int_from_cons(c);
        
        // calculate the proper momentum to add
        const double mu = get_mean_molecular_weight(c->prim[ZZZ]);
        const double n = c->prim[RHO] / (mu * m_proton);

        const double metallicity_solar = 0.02;
        double f_Z;
        if ((c->prim[ZZZ]/metallicity_solar) > .01)
        {
            f_Z = std::pow(c->prim[ZZZ]/metallicity_solar, -.14);
        }
        else
        {
            f_Z = 2;
        }


        const double p_t = 4.8e5
            * std::pow(E_blast/1e51, 13./14.)
            * std::pow(n, -1./7.)
            * std::pow(f_Z, 3./2.)
            * 1e5 // km/s in terms of cm/s
            * M_sun;
        printf("p_t           (cell %d) = %.3e M_sun km / s\n",
            i + n_guard_cell,
            p_t / (1e5* M_sun)
            );

        const double dp_prime_prime = dp_prime 
            * std::min(
                std::sqrt(1 + (c->cons[DDD]/(M_blast*weight))),
                p_t / p_blast
            );
        printf("p_prime_prime (cell %d) = %.3e M_sun km / s\n",
            i + n_guard_cell,
            dp_prime_prime / (1e5* M_sun)
            );



        // Apply ejected conservatives

        c->cons[DDD]    += dM;
        cons_added[DDD] += dM;

        c->cons[SRR]    += dp_prime_prime;
        cons_added[SRR] += dp_prime_prime;

        c->cons[TAU]    += dE;
        cons_added[TAU] += dE;

        c->cons[ZZZ]    += dMZ;
        cons_added[ZZZ] += dMZ;

        // possible adjust for cooling beyond cooling radius

        const double R_cool = 28.4 * pc_unit
            * std::pow(n, -3./7.)
            * std::pow(E_blast / 1e51, 2./7.)
            * f_Z;
        if (c->riph > R_cool)
        {
            printf("beyond R_cool for i=%d\n", i+n_guard_cell);
            printf("r      = %f pc\n", c->riph / pc_unit);
            printf("R_cool = %f pc\n", R_cool / pc_unit);
            printf("n = %e cm^-3 \n", n);
            printf("Z = %f\n", c->prim[ZZZ]);
            printf("f_Z = %e\n", f_Z);

            const double E_int_new_tmp = E_int_from_cons(c);

            const double dE_int = E_int_new_tmp - E_int_old;
            assert(dE_int > 0);

            const double dE_int_corrected = dE_int 
                * std::pow(c->riph / R_cool, -6.5);

            c->cons[TAU] += dE_int_corrected - dE_int;
        }


    }

    // correct for round-off error residuals
    struct cell * c_inner = &(theDomain->theCells[n_guard_cell]);
    c_inner->cons[DDD] += M_blast   - cons_added[DDD];
    c_inner->cons[TAU] += E_blast   - cons_added[TAU];
    c_inner->cons[ZZZ] += M_blast_Z - cons_added[ZZZ];

    // post conditions
    if (std::abs(1 - (cons_added[DDD] / M_blast))   > .05 )
    {
        printf("M_blast         = %f M_sun\n", M_blast / M_sun);
        printf("cons_added[DDD] = %f M_sun\n", cons_added[DDD] / M_sun);

        assert(std::abs(1 - (cons_added[DDD] / M_blast))   < .05 );
    }

    if (std::abs(1 - (cons_added[TAU] / E_blast))   > .05 )
    {
        printf("E_blast         = %e\n", E_blast);
        printf("cons_added[TAU] = %e\n", cons_added[TAU]);

        assert(std::abs(1 - (cons_added[TAU] / E_blast))   < .05 );
    }

    if (std::abs(1 - (cons_added[ZZZ] / M_blast_Z))   > .05 )
    {
        printf("M_blast_Z       = %f M_sun\n", M_blast_Z / M_sun);
        printf("cons_added[ZZZ] = %f M_sun\n", cons_added[ZZZ] / M_sun);

        assert(std::abs(1 - (cons_added[ZZZ] / M_blast_Z))   < .05 );
    }
    assert(std::abs(1 - (cons_added[TAU] / E_blast))   < .05 );
    assert(std::abs(1 - (cons_added[ZZZ] / M_blast_Z)) < .05 );
    std::cout.flush();

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

