
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
#include "Output/ascii.H" // count_lines_in_file

#include "blast.H"

using namespace boost;
using namespace boost::filesystem;

void add_winds( struct domain * theDomain , const double dt)
{
    // ============================================= //
    //
    //  Add a wind to account for pre-supernova mass loss
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
    //
    // ============================================= //

    add_winds_constant( theDomain, dt );

    // now we need to background within substep()
    // since we changed the cons, but not the prims
    calc_prim( theDomain );
    boundary( theDomain );


}

void add_winds_constant( struct domain * theDomain , const double dt ,
                         const double wind_velocity )
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
    //     - wind_velocity  - injection velocity
    //                      -  default: 1000 km/s
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


    const int n_guard_cell = 1;
    struct cell * blast_cell = &(theDomain->theCells[n_guard_cell]);


    for ( unsigned int i=0 ; i <theDomain->SNe.size() ; ++i )
    {
        supernova SN = theDomain->SNe[i];
        const double mass_loss = (dt / SN.lifetime) 
                                * (SN.mass - SN.mass_ejecta);

        blast_cell->cons[DDD] += mass_loss;
        blast_cell->cons[SRR] += mass_loss * wind_velocity;
        blast_cell->cons[TAU] += mass_loss * .5 * std::pow(wind_velocity, 2); 
        blast_cell->cons[ZZZ] += mass_loss * theDomain->metallicity;

    }


    return;

}

double get_ejecta_mass( const double M_initial )
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
                                            11.837161, 12.628000, 13.283933,
                                            13.973949, 14.867287, 15.229942,
                                            15.757695, 15.604063, 16.343219,
                                            16.144743, 15.389504, 14.995688,
                                            15.197527, 15.178074, 14.443274,
                                            14.193992, 14.117690, 13.952423,
                                            13.509647, 13.313397, 13.354588,
                                            15.109093, 10.222817,  8.271961,
                                             5.883033,  4.969768,  4.579023,
                                             4.227612,  4.386078 };
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

double get_ejecta_mass_Z( const double M_initial )
{
    // ============================================= //
    //
    //  Given an initial mass, find the ejecta mass (mass)
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

    c->cons[DDD] += M_blast;
    c->cons[TAU] += E_blast;

    // for now assume that the metallicity of the ejecta is the same
    // as the background metallicity
    // later we'll want to actually inject metals
    c->cons[ZZZ] += M_blast_Z;


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


std::vector<supernova> get_SNe( const double cluster_mass,
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
    //     - metallicity    - metallicity mass fraction (e.g. solar = 0.02)
    //                      - Default = 0.02
    //     - seed           - seeds the stochasticity;
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
            tmp.mass_ejecta   = get_ejecta_mass(  mass) * M_sun; // [g]
            tmp.mass_ejecta_Z = get_ejecta_mass_Z(mass) * M_sun; // [g]
            SNe.push_back(tmp);
        }
    }

    // sort in reverse order, so we can pop off SNe as they happen
    std::sort(SNe.rbegin(), SNe.rend(), sort_by_lifetime);

    std::cout << "Num SNe: " << SNe.size() << std::endl;
    std::cout << "Cluster mass: " << cluster_mass << std::endl;

    return SNe;
}

std::vector<supernova> read_SNe( const std::string filename)
{

    std::vector<supernova> SNe;
    supernova SN_tmp;

    int nL = count_lines_in_file(filename) - 1;
    if ( nL < 0 ) return SNe;

    double SN_time;
    double SN_mass;
    double SN_mass_ejecta;
    double SN_mass_ejecta_Z;


    FILE * pFile = fopen(filename.c_str(),"r");
    char tmp[1024];
    fgets(tmp, sizeof(tmp), pFile); // header line

    for( int l=0 ; l<nL ; ++l )
    {
        fscanf(pFile,"%le %le %le %le\n",
                &SN_time, &SN_mass, &SN_mass_ejecta, &SN_mass_ejecta_Z);

        SN_tmp.mass             = SN_mass;
        SN_tmp.mass_ejecta      = SN_mass_ejecta;
        SN_tmp.mass_ejecta_Z    = SN_mass_ejecta_Z;
        SN_tmp.lifetime         = SN_time;

        SNe.push_back(SN_tmp);
    }
    fclose(pFile);

    std::sort(SNe.rbegin(), SNe.rend(), sort_by_lifetime);

    return SNe;
}

bool sort_by_lifetime( const supernova &a, const supernova &b)
{
    return (a.lifetime < b.lifetime);
}

