
#include <assert.h>
#include <cmath>

#include "structure.H"

#include "boundary.H" // boundary
#include "cooling.H"
#include "misc.H" // calc_dr, calc_prim, radial_flux, add_source, etc
#include "timestep.H"
#include "Initial/initial_conditions.H"
#include "Hydro/euler.H" // E_int_from_cons






void update_initial_energies( struct domain * theDomain )
{
    // ============================================= //
    //
    //  Only used for Keller subgrid conduction method
    //  Updates a tracker of the initial kinetic and thermal energies
    //    for each cell.
    //  Also updates the energy split between the hot and cold subgrid components
    //  Also updates the mass split between hot and cold
    //  
    //
    //  Inputs:
    //     - theDomain    - the standard domain struct used throughout
    //
    // ============================================= //
    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;   

    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        struct cell * c = &(theCells[i]);

        // ----------- Pre-conditions ------------- //
        if (c->multiphase)
        {
            verify_multiphase_conditions(c, "update_initial_energies()", "pre");

            if (E_int_from_cons(c->cons_hot) < 0)
            {
                printf("------ERROR in update_initial_energies() pre-conditions------- \n");
                printf("hot gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_hot) = %e \n", E_int_from_cons(c->cons_hot));
                assert( E_int_from_cons(c->cons_hot) > 0 );
            }

            if (E_int_from_cons(c->cons_cold) < 0)
            {
                printf("------ERROR in update_initial_energies() pre-conditions------- \n");
                printf("cold gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_cold) = %e \n", E_int_from_cons(c->cons_cold));
                assert( E_int_from_cons(c->cons_cold) > 0 );
            }
        }

        if(c->multiphase)
        {
            c->E_int_initial = E_int_from_cons(c->cons);
            c->E_kin_initial = E_kin_from_cons(c->cons);

            c->x_hot  = c->cons_hot[DDD]  / c->cons[DDD];
            c->x_cold = c->cons_cold[DDD] / c->cons[DDD];

            c->y_hot  = E_int_from_cons(c->cons_hot)  / E_int_from_cons(c->cons);
            c->y_cold = E_int_from_cons(c->cons_cold) / E_int_from_cons(c->cons);

            c->z_hot  = c->cons_hot[ZZZ]  / c->cons_hot[DDD];
            c->z_cold = c->cons_cold[ZZZ] / c->cons_cold[DDD];
        }

        // ----------- Post-conditions ------------- //
        if (c->multiphase)
        {
            const double rel_tol = 1e-5; // relative tolerance for float comparisons

            if ((c->x_hot < 0) || (c->x_hot > 1+rel_tol))
            {
                printf("------ERROR in update_initial_energies() post-conditions ------- \n");
                printf("hot gas mass fraction not in [0, 1+rel_tol]\n");
                printf("c->x_hot          = %e \n", c->x_hot);
                printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
                printf("c->cons[DDD]      = %e \n", c->cons[DDD]);
                assert( (c->x_hot > 0) && (c->x_hot < 1+rel_tol) );
            }

            if ((c->x_cold < 0) || (c->x_cold > 1+rel_tol))
            {
                printf("------ERROR in update_initial_energies() post-conditions ------- \n");
                printf("cold gas mass fraction not in [0, 1+rel_tol]\n");
                printf("c->x_cold         = %e \n", c->x_cold);
                printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
                printf("c->cons[DDD]      = %e \n", c->cons[DDD]);
                assert( (c->x_cold > 0) && (c->x_cold < 1+rel_tol) );
            }

            if ((c->y_hot < 0) || (c->y_hot > 1+rel_tol))
            {
                printf("------ERROR in update_initial_energies() post-conditions ------- \n");
                printf("hot gas energy fraction not in [0, 1+rel_tol]\n");
                printf("c->y_hot                      = %e \n", c->y_hot);
                printf("E_int_from_cons(c->cons_hot)  = %e \n", E_int_from_cons(c->cons_hot));
                printf("E_int_from_cons(c->cons_cold) = %e \n", E_int_from_cons(c->cons_cold));
                printf("E_int_from_cons(c->cons)      = %e \n", E_int_from_cons(c->cons));
                printf("c->cons_hot[TAU]              = %e \n", c->cons_hot[TAU]);
                printf("c->cons[TAU]                  = %e \n", c->cons[TAU]);
                assert( (c->y_hot > 0) && (c->y_hot < 1+rel_tol) );
            }

            if ((c->y_cold < 0) || (c->y_cold > 1+rel_tol))
            {
                printf("------ERROR in update_initial_energies() post-conditions ------- \n");
                printf("cold gas energy fraction not in [0, 1+rel_tol]\n");
                printf("c->y_cold                      = %e \n", c->y_cold);
                printf("E_int_from_cons(c->cons_hot)   = %e \n", E_int_from_cons(c->cons_hot));
                printf("E_int_from_cons(c->cons_cold)  = %e \n", E_int_from_cons(c->cons_cold));
                printf("E_int_from_cons(c->cons)       = %e \n", E_int_from_cons(c->cons));
                printf("c->cons_cold[TAU]              = %e \n", c->cons_cold[TAU]);
                printf("c->cons[TAU]                   = %e \n", c->cons[TAU]);
                assert( (c->y_cold > 0) && (c->y_cold < 1+rel_tol) );
            }
        }
    }
}

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
    //     - last_step    - is this the last substep within an overall step?
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

    update_initial_energies( theDomain );

    EnergyChecker energy_checker_flux(theDomain, "radial_flux");
    radial_flux( theDomain , dt );
    energy_checker_flux.check( theDomain );

    EnergyChecker energy_checker_source(theDomain, "add_source");
    add_source( theDomain , dt , cooling );
    if(!theDomain->theParList.With_Cooling) energy_checker_source.check( theDomain );

    EnergyChecker energy_checker_move_cells(theDomain, "move_cells");
    if( first_step ) move_cells( theDomain , dt );
    calc_dr( theDomain );
    energy_checker_move_cells.check( theDomain );

    // fix_negative_energies( theDomain );

    calc_prim( theDomain );

    if( last_step )
    {
        EnergyChecker energy_checker_amr(theDomain, "AMR");
        AMR( theDomain );
        energy_checker_amr.check( theDomain );


        theDomain->R_shock = ICs->find_shock( theDomain );
        ICs->possibly_extend_grid( theDomain );
        // check_multiphase( theDomain );

    }
    EnergyChecker energy_checker_boundary(theDomain, "boundary");
    boundary( theDomain );
    energy_checker_boundary.check( theDomain );


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
        if ( c->multiphase )
        {
            memcpy( c->RKcons_hot  , c->cons_hot  , NUM_Q*sizeof(double) );
            memcpy( c->RKcons_cold , c->cons_cold , NUM_Q*sizeof(double) );

        }
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


