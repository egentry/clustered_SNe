
#include "structure.H"

#include "boundary.H" // boundary
#include "cooling.H"
#include "misc.H" // calc_dr, calc_prim, radial_flux, add_source, etc
#include "timestep.H"
#include "Initial/initial_conditions.H"


void substep( struct domain * theDomain , double RK , 
              double dt , int first_step , int last_step ,
              Initial_Conditions * ICs ,
              Cooling * cooling )
{

    // ============================================= //
    //
    //  Evolve for a single SUB-STEP 
    //
    //  Inputs:
    //     - theDomain    - the standard domain struct used throughout
    //     - RK           - Runge-Kutta coefficient
    //                    - see note within adjustRK
    //     - dt           - time elapsed within substep
    //     - first_step   - is this the first substep within an overall step?
    //     - last_step    - is this the 
    //
    //  Outputs:
    //       None
    //
    //  Side effects:
    //     - Every cell:
    //       - Every subset - overwrites cons, prims
    //                      - can overwrite dr, but cells only move during first_step
    //       - Every first_step - overwrites dr
    //     - Boundary cells:
    //          - Can overwrite whatever happens within boundary()
    //
    //  Notes:
    //    - Does NOT change wiph, except if in boundary()
    //       - Might only be 1st order accurate in time
    //         if wiph and riph not integrated in a RK fashion?
    //
    // ============================================= //

    adjust_RK_cons( theDomain , RK );

    radial_flux( theDomain , dt );
    add_source( theDomain , dt , cooling );

    if( first_step ) move_cells( theDomain , dt );
    calc_dr( theDomain );

    fix_negative_energies( theDomain );

    calc_prim( theDomain );

    if( last_step )
    {
        AMR( theDomain );
        ICs->possibly_extend_grid( theDomain );

    }
    boundary( theDomain );

}


void timestep( struct domain * theDomain , const double dt,
              Initial_Conditions * ICs ,
              Cooling * cooling )
{
    // ============================================= //
    //
    //  Evolve a full step (possibly with Runge-Kutta substeps)
    //
    //  Inputs:
    //     - theDomain    - the standard domain struct used throughout
    //     - dt           - time elapsed within substep
    //
    //  Outputs:
    //       None
    //
    //  Side effects:
    //     - Overwrites RKcons for each cell
    //     - Also has all the side effects in substep()
    //
    //  Notes:
    //
    // ============================================= //

    struct cell * theCells = theDomain->theCells;
    int Nr = theDomain->Nr;

    int i;

    for( i=0 ; i<Nr ; ++i )
    {
        struct cell * c = &(theCells[i]);
        memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
    }

    // for choosing proper Runge-Kutta coefficients for higher order schemes,
    //   see: Gottlieb & Chi-Wang Shu (1998)
    //        "Total Variation Diminishing Runge-Kutta Schemes"
    //        Mathematics of Computation
    //        http://www.ams.org/journals/mcom/1998-67-221/S0025-5718-98-00913-2/S0025-5718-98-00913-2.pdf
    if ( theDomain->theParList.RK2 )
    {
        substep( theDomain , 0.0 ,     dt , 1 , 0 , ICs , cooling );
        substep( theDomain , 0.5 , 0.5*dt , 0 , 1 , ICs , cooling );      
    }
    else
    {
        substep( theDomain , 0.0 ,     dt , 1 , 1 , ICs , cooling );
    }

    theDomain->t += dt;   
    theDomain->count_steps += 1;

}


