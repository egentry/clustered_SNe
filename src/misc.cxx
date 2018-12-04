
#include <string.h>
#include <assert.h>
#include <cmath> // std::abs, std::isfinite

extern "C" {
#include <grackle.h>
}

#include "geometry.H" // get_dA, get_dV
#include "structure.H"
#include "constants.H"
#include "misc.H" 
#include "plm.H"
#include "riemann.H"
#include "Hydro/euler.H" // cons2prim, mindt, E_int_from_*, calc_T



double getmindt( const struct domain * theDomain )
{

    // ============================================= //
    //
    //  Calculates the overall timestep
    //
    //  Inputs:
    //     - theDomain    - the standard domain struct used throughout
    //
    //  Returns:
    //    - dt            - the overall timestep to be used
    //
    //  Side effects:
    //    None 
    //
    //  Notes:
    //
    // ============================================= //

    const struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;

    double dt = 1e100;
    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        const int im = i-1;
        const struct cell * c = &(theCells[i]);
        const double dr = c->dr;
        const double r = c->riph-.5*dr;
        const double wm = theCells[im].wiph;
        const double wp = c->wiph;
        const double w = .5*(wm+wp);

        double dt_temp = mindt( c->prim , w , r , dr );

        if( dt > dt_temp ) dt = dt_temp;
    }
    dt *= theDomain->theParList.CFL; 

    return dt;
}

double get_mean_molecular_weight(const double Z)
{
    // assumes atomic, fully ionized species
    const double Y = .23; // helium fraction
    const double X = 1 - Y - Z; // hydrogen mass fraction

    const double mu = 1. / (2*X + .75*Y + .5*Z); // mean molecular weight

    return mu;
}

void set_wcell( struct domain * theDomain )
{

    // ============================================= //
    //
    //  Calculates and sets velocities of cell boundaries
    //
    //  Inputs:
    //     - theDomain    - the standard domain struct used throughout
    //
    //  Returns:
    //    void
    //
    //  Side effects:
    //     - overwrites the wiph variable for each cell
    //
    //  Notes:
    //    - Only updated once per timestep,
    //      even if we're using higher order Runge Kutta sub steps
    //    - "w" is defined as the velocity as the *outer* boundary of cell
    //       - w is not set for outermost cell
    //
    // ============================================= //

    struct cell * theCells = theDomain->theCells;
    const int mesh_motion = theDomain->theParList.Mesh_Motion;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;

    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        struct cell * cL = &(theCells[i]);  
        double w = 0.0;
        if( mesh_motion )
        {
            // if Lagrangian, average the cell-centered velocities to find the edge velocity in between
            const struct cell * cR = &(theCells[i+1]);
            const double wL = cL->prim[VRR];
            const double wR = cR->prim[VRR];
            w = .5*(wL + wR); 
            // why this special case for the outside edge of the innermost cell?
            // if( i==0 ) w = wR*(cR->riph - .5*cR->dr)/(cL->riph);//0.0;//2./3.*wR;
            // if(i==0) w=0;
        }
        cL->wiph = w;
    }
    // Boundary condition: innermost cell doesn't move
    theCells[0].wiph = 0;

}

void adjust_RK_cons( struct domain * theDomain , const double RK )
{

    // ============================================= //
    //
    //  Sets the cons variable for a given substep,
    //    using a combination of cons and RKcons
    //       -   cons: current conservative variables after previous subset
    //       - RKcons: conservative variables at start of this overall step
    //
    //  Inputs:
    //     - theDomain    - the standard domain struct used throughout
    //     - RK           - the fraction of cons that should be from the 
    //                      original values at start of this substep
    //                      (i.e. RK=1.0 overwrites current cons)
    //
    //  Returns:
    //    void
    //
    //  Side effects:
    //     - overwrites the cons variable
    //
    //  Notes:
    //     - Assumes RKcons was properly set to valid cons values
    //       at the start of this timestep (not this sub-step)
    //
    //     - for choosing proper Runge-Kutta coefficients for higher order schemes,
    //       see: Gottlieb & Chi-Wang Shu (1998)
    //            "Total Variation Diminishing Runge-Kutta Schemes"
    //            Mathematics of Computation
    //            http://www.ams.org/journals/mcom/1998-67-221/S0025-5718-98-00913-2/S0025-5718-98-00913-2.pdf
    //       for 2nd order, see equation 2.4
    //
    // ============================================= //

    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;

    for( int i=0 ; i<Nr ; ++i )
    {
        struct cell * c = &(theCells[i]);


        // ----------- Pre-conditions ------------- //

        if (c->multiphase)
        {
            verify_multiphase_conditions(c, "adjust_RK_cons()", "pre");

            if (E_int_from_cons(c->cons_hot) < 0)
            {
                printf("------ERROR in adjust_RK_cons() pre-conditions------- \n");
                printf("hot gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_hot) = %e \n", E_int_from_cons(c->cons_hot));
                printf("i = %d\n", i);
                assert( E_int_from_cons(c->cons_hot) > 0 );
            }

            if (E_int_from_cons(c->cons_cold) < 0)
            {
                printf("------ERROR in adjust_RK_cons() pre-conditions------- \n");
                printf("cold gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_cold) = %e \n", E_int_from_cons(c->cons_cold));
                printf("i = %d\n", i);
                assert( E_int_from_cons(c->cons_cold) > 0 );
            }
        }

        for( int q=0 ; q<NUM_Q ; ++q )
        {
            c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];

            if (c->multiphase)
            {
                c->cons_hot[q]  = (1.-RK)*c->cons_hot[q]  + RK*c->RKcons_hot[q];
                c->cons_cold[q] = (1.-RK)*c->cons_cold[q] + RK*c->RKcons_cold[q];
            }
        }

        // ======== Verify post-conditions ========= //
        #ifndef NDEBUG
        double prim_tmp[NUM_Q];
        const double rp = c->riph;
        const double rm = rp-c->dr;
        const double dV = get_dV( rp , rm );
        cons2prim( c->cons , prim_tmp , dV );

        if(prim_tmp[PPP] < theDomain->theParList.Pressure_Floor)
        {
            printf("------ ERROR in adjust_RK_cons()------- \n");
            printf("pressure should be above pressure floor! \n");
            printf("pressure       = %e in cell %i \n", prim_tmp[PPP], i);
            printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
            printf("rp = %e \n", rp);
            printf("rm = %e \n", rm);
            printf("dr = %e \n", c->dr);
            assert(prim_tmp[PPP] < theDomain->theParList.Pressure_Floor);
        }
        for( int q=0 ; q<NUM_Q ; ++q)
        {
            if( !std::isfinite(prim_tmp[q]) )
            {
               printf("------ ERROR in adjust_RK_cons()------- \n");
               printf("prim[%d] = %e in cell %d \n", q, prim_tmp[q], i);
               printf("rp = %e \n", rp);
               printf("rm = %e \n", rm);
               printf("dr = %e \n", c->dr);
               assert( std::isfinite(prim_tmp[q]) );
            }
        }


        if (c->multiphase)
        {
            verify_multiphase_conditions(c, "adjust_RK_cons()", "post");

            if (E_int_from_cons(c->cons_hot) < 0)
            {
                printf("------ERROR in adjust_RK_cons() post-conditions------- \n");
                printf("hot gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_hot) = %e \n", E_int_from_cons(c->cons_hot));
                printf("i = %d\n", i);
                assert( E_int_from_cons(c->cons_hot) > 0 );
            }

            if (E_int_from_cons(c->cons_cold) < 0)
            {
                printf("------ERROR in adjust_RK_cons() post-conditions------- \n");
                printf("cold gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_cold) = %e \n", E_int_from_cons(c->cons_cold));
                printf("i = %d\n", i);
                assert( E_int_from_cons(c->cons_cold) > 0 );
            }
        }

        #endif
    }
}

void move_cells( struct domain * theDomain , const double dt)
{

    // ============================================= //
    //
    //  Moves cell boundaries
    //
    //  Inputs:
    //    - theDomain    - the standard domain struct used throughout
    //    - RK           - not used
    //
    //  Returns:
    //    void
    //
    //  Side effects:
    //    - overwrites the riph variable for each cell
    //
    //  Notes:
    //    - Only updated once per timestep (on "first_step")
    //    - riph is defined as the velocity as the *outer* boundary of cell
    //    - Assumed that the boundary cell w's are properly set
    //      in set_wcell
    //
    // ============================================= //

    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    for( int i=0 ; i<Nr ; ++i )
    {
        struct cell * c = &(theCells[i]);

        // ======== Verify pre-conditions ========= //
        #ifndef NDEBUG
         
         double prim_tmp[NUM_Q];
         double rp = c->riph;
         double rm = rp-c->dr;
         double dV = get_dV( rp , rm );
         cons2prim( c->cons , prim_tmp , dV );
         if(prim_tmp[PPP] < theDomain->theParList.Pressure_Floor)
         {
            printf("------ ERROR in move_cells() before moving------- \n");
            printf("pressure should be above pressure floor! \n");
            printf("pressure       = %e in cell %i \n", prim_tmp[PPP], i);
            printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
            printf("rp = %e \n", c->riph);
            printf("rm = %e \n", c->riph - c->dr);
            printf("dr = %e \n", c->dr);
            assert( prim_tmp[PPP] < theDomain->theParList.Pressure_Floor );
         }
        #endif

        // ======== Move Cells ========= //

        c->riph += c->wiph*dt;

        // ======== Verify post-conditions ========= //
        #ifndef NDEBUG
        rp = c->riph;
        rm = rp-c->dr;
        dV = get_dV( rp , rm );
        cons2prim( c->cons , prim_tmp , dV );
        if(prim_tmp[PPP] < theDomain->theParList.Pressure_Floor)
        {
            printf("------ ERROR in move_cells() after moving------- \n");
            printf("pressure should be above pressure floor! \n");
            printf("pressure       = %e in cell %i \n", prim_tmp[PPP], i);
            printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
            printf("rp = %e \n", rp);
            printf("rm = %e \n", rm);
            printf("dr = %e \n", c->dr);
            assert(prim_tmp[PPP] < theDomain->theParList.Pressure_Floor);
        }
        for( int q=0 ; q<NUM_Q ; ++q)
        {
            if( !std::isfinite(prim_tmp[q]) )
            {
                printf("------ ERROR in move_cells()------- \n");
                printf("non-finite prim[q] \n");
                printf("prim[%d] = %e in cell %d \n", q, prim_tmp[q], i);
                printf("rp = %e \n", rp);
                printf("rm = %e \n", rm);
                printf("dr = %e \n", c->dr);
                assert( std::isfinite(prim_tmp[q]) );
            }
        }
        #endif
    }
}

void calc_dr( struct domain * theDomain )
{

    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;

    for( int i=1 ; i<Nr ; ++i )
    {
        int im = i-1;
        const double rm = theCells[im].riph;
        const double rp = theCells[i ].riph;
        const double dr = rp-rm;
        theCells[i].dr = dr;
        assert(dr > 0);

    }
    // Boundary condition: innermost cell extends to r=0
    theCells[0].dr = theCells[0].riph;

}

void fix_negative_energies( struct domain * theDomain )
{
    // Only TOTAL energy is tracked in our riemann solver
    //    - this can lead to a problem if kinetic energy exceeds total energy
    //    - in order to compensate, internal energy (and pressure) will become negative
    // To fix this problem:
    //    - we specifically evolve internal energy using adiabatic assumptions
    //    - if the adiabatic internal energy differs from the numeric internal energy by too much
    //       then we use the adiabatic energy rather than the numeric energy

    // 

    const double tolerance = .9; // relative tolerance allowed on energy error before switching to adiabatic

    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;
    const double gamma = theDomain->theParList.Adiabatic_Index;

    for( int i=Ng ; i<Nr-Ng ; ++i ) // don't extend to outermost zone; that already doesn't conserve energy adiabatically
    {  
        struct cell * c = &(theCells[i]);

        const double rp = c->riph;
        const double rm = rp - c->dr;
        const double dV = get_dV( rp , rm );

        const double E_int_adiabatic = c->E_int_old * pow(dV / c->dV_old, 1-gamma);
        
        double dE = c->dE_cool;
        if(theDomain->theParList.with_turbulent_diffusion)
        {
            dE += c->dE_diffusion;
        }
        const double E_int_approx = E_int_adiabatic + dE;


        const double E_int_numeric = E_int_from_cons( c->cons );

        if ( (((E_int_approx - E_int_numeric) / std::abs(E_int_approx)) > tolerance) 
                || (E_int_numeric < 0) ) 
        {

            const double mass  = c->cons[DDD];
            const double vr    = c->cons[SRR] / mass;
            // overwrite the energy, so that internal energy is strictly positive
            printf("------ Applying energy fix for cell %d  ------- \n", i);
            printf("time            = %e \n", theDomain->t);
            printf("velocity        = %e \n", vr);
            printf("E_int_old       = %e \n", c->E_int_old);
            printf("dV_old          = %e \n", c->dV_old);
            printf("dV              = %e \n", dV);
            printf("P_approx        = %e \n", E_int_approx * (gamma-1)/dV);
            printf("E_int_adiabatic = %e \n", E_int_adiabatic);
            printf("E_int_approx    = %e \n", E_int_approx);
            printf("dE_cool         = %e \n", c->dE_cool);
            printf("E_int_numeric   = %e \n", E_int_numeric);

            if ( E_int_approx > 0 )
            {
                c->E_int_old = E_int_approx;
            }
            else
            {
                c->E_int_old = E_int_adiabatic;
            }
            c->cons[TAU] = c->E_int_old + .5*mass*vr*vr;
            if(c->multiphase)
            {
                // break, since I haven't thought about how to implement 
                // dual energy formalism for multiphase subgrid
                assert(0);
            }
        } 
        else
        {
            c->E_int_old  = E_int_numeric;
        }
        
        c->dV_old = dV;


    }

   return;
}


void calc_prim( struct domain * theDomain )
{

    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;

    for( int i=0 ; i<Nr ; ++i )
    {
        struct cell * c = &(theCells[i]);
        const double rp = c->riph;
        const double rm = rp-c->dr;
        const double dV = get_dV( rp , rm );
        cons2prim( c->cons , c->prim , dV );
        if (c->multiphase)
        {
            calc_multiphase_prim( c, c->prim_hot , c->prim_cold ,
                                     &(c->V_hot) , &(c->V_cold) );   
        }


        // ======== Verify post-conditions ========= //
        #ifndef NDEBUG
        if(c->prim[PPP] < theDomain->theParList.Pressure_Floor)
        {
            printf("------ ERROR in calc_prim()------- \n");
            printf("pressure should be above pressure floor! \n");
            printf("pressure       = %e in cell %i \n", c->prim[PPP], i);
            printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
            printf("rp = %e \n", rp);
            printf("rm = %e \n", rm);
            printf("dr = %e \n", c->dr);
            assert(c->prim[PPP] < theDomain->theParList.Pressure_Floor);
        }
        for( int q=0 ; q<NUM_Q ; ++q)
        {
            if( !std::isfinite(c->prim[q]) )
            {
                printf("------ ERROR in calc_prim()------- \n");
                printf("prim[%d] = %e in cell %d \n", q, c->prim[q], i);
                printf("rp = %e \n", rp);
                printf("rm = %e \n", rm);
                printf("dr = %e \n", c->dr);
                assert( std::isfinite(c->prim[q]) );
            }
        }
        #endif
    }
}


void radial_flux( struct domain * theDomain , const double dt )
{

    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;

    plm( theDomain );
    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        struct cell * cL = &(theCells[i]);
        struct cell * cR = &(theCells[i+1]);
        const double r = cL->riph;
        const double dA = get_dA(r); 
        if(cL->multiphase || cR->multiphase)
        {
            printf("i = %d (radial_flux() -- multiphase))\n", i);
        }
        riemann( cL , cR , r , dA , dt );
        // artificial_conduction( cL , cR , dA , dt ); // Noh's artificial conduction
        // subgrid_thermal_conduction(cL, dt); // Keller's inter-cell conduction
    }

    if(theDomain->theParList.with_physical_conduction)
    {
        thermal_conduction_implicit( theDomain , dt );
        // thermal_conduction_explicit( theDomain , dt );
    }

    if(theDomain->theParList.with_turbulent_diffusion)
    {
        // turbulent_diffusion_implicit_mass_density(  theDomain , dt );
        // turbulent_diffusion_implicit_momentum_density(  theDomain , dt );
        // turbulent_diffusion_implicit_E_tot_density( theDomain , dt );
        // turbulent_diffusion_implicit_metal_density( theDomain , dt );

        EnergyChecker energy_checker_turb(theDomain, "turbulent_diffusion_implicit");
        turbulent_diffusion_implicit( theDomain , dt );
        energy_checker_turb.check( theDomain );
    }

}


void add_source( struct domain * theDomain , const double dt , 
                 Cooling * cooling )
{

    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;

    double net_energy_from_cooling = 0;

    for( int i=Ng ; i<Nr ; ++i )
    {
        struct cell * c = &(theCells[i]);
        const double rp = c->riph;
        double rm = rp-c->dr;
        const double dV = get_dV(rp,rm);
        if (i==1)
        {
            rm = 0; // boundary condition -- don't change dV to match this rm
        }
        if(c->multiphase)
        {
            printf("i = %d (in source() -- only for multiphase)\n", i);
        }
        source( c , rp , rm , dV , dt , theDomain->R_shock , cooling );
        net_energy_from_cooling += c->dE_cool;
        if((net_energy_from_cooling>0) & (c->dE_cool > 0))
        {
            printf("net heating about threshold -- cell %d; time %e [Myr]; dt = %e \n", 
                i, theDomain->t / (1e6 * yr), dt);
            printf("dE_cool/dt = %e [erg / Myr]\n", c->dE_cool / dt * (1e6 * yr));
            printf("total energy added by heating = %e \n", theDomain->energy_added_by_cooling);
            printf("net so far = %e \n", net_energy_from_cooling);
            printf("dE_cool = %e \n", c->dE_cool);
            printf("E(before cool) = %e \n", c->cons[TAU] - c->dE_cool);
            printf("Z = %e \n", c->prim[ZZZ]);
            printf("density = %e \n", c->prim[RHO]);
            printf("velocity = %e \n", c->prim[VRR]);
            printf("pressure = %e  \n", c->prim[VRR]);

            fflush(stdout);
        }

        if((c->dE_cool>0) & ((c->dE_cool / (c->cons[TAU] - c->dE_cool)) > .01))
        {
            printf("significant fraction of net heating -- cell %d \n", i);
            printf("dE_cool        = %e \n", c->dE_cool);
            printf("E(before cool) = %e \n", c->cons[TAU] - c->dE_cool);
            printf("dE = %e %%\n", c->dE_cool / (c->cons[TAU] - c->dE_cool) *100);
            fflush(stdout);
        }


        // ======== Verify post-conditions ========= //
        #ifndef NDEBUG
        if(c->prim[PPP] < theDomain->theParList.Pressure_Floor)
        {
            printf("------ ERROR in add_source()------- \n");
            printf("pressure should be above pressure floor! \n");
            printf("pressure       = %e in cell %i \n", c->prim[PPP], i);
            printf("pressure floor = %e \n", theDomain->theParList.Pressure_Floor);
            printf("rp = %e \n", rp);
            printf("rm = %e \n", rm);
            printf("dr = %e \n", c->dr);
            assert(c->prim[PPP] < theDomain->theParList.Pressure_Floor);
        }
        for( int q=0 ; q<NUM_Q ; ++q)
        {
            if( !std::isfinite(c->prim[q]) )
            {
                printf("------ ERROR in add_source()------- \n");
                printf("prim[%d] = %e in cell %d \n", q, c->prim[q], i);
                printf("rp = %e \n", rp);
                printf("rm = %e \n", rm);
                printf("dr = %e \n", c->dr);
                assert( std::isfinite(c->prim[q]) );
            }
        }
        #endif
    }
    if(net_energy_from_cooling > 0)   
    {
        theDomain->energy_added_by_cooling += net_energy_from_cooling;
    }

}


void longandshort( const struct domain * theDomain , 
                   double * L , double * S , 
                   int * iL , int * iS )
{ 

    const struct cell * theCells = theDomain->theCells;
    const int    Nr    = theDomain->Nr;
    const int    Ng    = theDomain->Ng;
    const double rmax  = theCells[Nr-1].riph;
    const double rmin  = theCells[0].riph;
    const int    Nr0   = theDomain->theParList.Num_R;
    const double dr0   = rmax / static_cast <double> (Nr0);
    const double dx0   = log(rmax/rmin)/Nr0;
    const int logscale = theDomain->theParList.LogZoning;

    double Long  = 0.0; 
    double Short = 0.0; 
    int iLong    = -1;
    int iShort   = -1;


    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        double l,s;
        const struct cell * c = &(theCells[i]);
        if( logscale )
        {
            const double dy = c->dr;
            const double dx = c->riph*dx0;
            l = dy/dx;
            s = dx/dy;
        }
        else
        {
            const double dy = c->dr;
            const double dx = dr0;
            l = dy/dx;
            s = dx/dy;
        }
        if( Long  < l ){ Long  = l; iLong  = i; } 
        if( Short < s ){ Short = s; iShort = i; } 
    }

    *iS = iShort;
    *iL = iLong;
    *S = Short;
    *L = Long;

}


void AMR( struct domain * theDomain )
{

    double L,S;
    int iL=0;
    int iS=0;
    longandshort( theDomain , &L , &S , &iL , &iS );

    // Thresholds for applying AMR
    const double MaxShort = theDomain->theParList.MaxShort;
    const double MaxLong  = theDomain->theParList.MaxLong;

    struct cell * theCells = theDomain->theCells;
    int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;

    if( S>MaxShort )
    {
        if (iS == 1)
        {
            // printf("Kill! iS = %d\n",iS);         
        }
        printf("KILL!  iS = %d\n",iS);

        int iSp = iS+1;
        int iSm = iS-1;
        int iMerge = iSp;
        if ( iS>Ng )
        {
            //Possibly shift iS backwards by 1 
            double drL = theCells[iSm].dr;
            double drR = theCells[iSp].dr;
            if( drL<drR )
            {
                --iS;
                --iSm;
                --iSp;
                iMerge -=2;
            }
        }

        printf("merge with i=%d\n", iMerge);


        struct cell * c  = &(theCells[iS]);
        const struct cell * cp = &(theCells[iSp]);

        if( (c->multiphase) || (cp->multiphase) )
        {
            printf("------ Before ------\n");
            printf("Left ------\n");

            printf("\n");

            printf("c->cons[DDD]       = %e \n", c->cons[DDD]);
            printf("c->cons[SRR]       = %e \n", c->cons[SRR]);
            printf("c->cons[TAU]       = %e \n", c->cons[TAU]);

            printf("\n");

        if(c->multiphase)
        {
            printf("c->cons_hot[DDD]   = %e \n", c->cons_hot[DDD]);
            printf("c->cons_hot[SRR]   = %e \n", c->cons_hot[SRR]);
            printf("c->cons_hot[TAU]   = %e \n", c->cons_hot[TAU]);

            printf("\n");

            printf("c->cons_cold[DDD]  = %e \n", c->cons_cold[DDD]);
            printf("c->cons_cold[SRR]  = %e \n", c->cons_cold[SRR]);
            printf("c->cons_cold[TAU]  = %e \n", c->cons_cold[TAU]);

            printf("\n");
        }

            printf("E_int (total)      = %e \n", E_int_from_cons(c->cons));
            if(c->multiphase)
            {
                printf("E_int (hot)        = %e \n", E_int_from_cons(c->cons_hot));
                printf("E_int (cold)       = %e \n", E_int_from_cons(c->cons_cold));
            }

            printf("\n");

            printf("v (total)          = %e \n", c->cons[SRR]      / c->cons[DDD]);
            if(c->multiphase)
            {
                printf("v (hot)            = %e \n", c->cons_hot[SRR]  / c->cons_hot[DDD]);
                printf("v (cold)           = %e \n", c->cons_cold[SRR] / c->cons_cold[DDD]);
            }

            printf("\n");

            printf("Right ------\n");
            printf("\n");

            printf("cp->cons[DDD]      = %e \n", cp->cons[DDD]);
            printf("cp->cons[SRR]      = %e \n", cp->cons[SRR]);
            printf("cp->cons[TAU]      = %e \n", cp->cons[TAU]);

            printf("\n");

            if(cp->multiphase)
            {
                printf("cp->cons_hot[DDD]  = %e \n", cp->cons_hot[DDD]);
                printf("cp->cons_hot[SRR]  = %e \n", cp->cons_hot[SRR]);
                printf("cp->cons_hot[TAU]  = %e \n", cp->cons_hot[TAU]);


                printf("\n");

                printf("cp->cons_cold[DDD] = %e \n", cp->cons_cold[DDD]);
                printf("cp->cons_cold[SRR] = %e \n", cp->cons_cold[SRR]);
                printf("cp->cons_cold[TAU] = %e \n", cp->cons_cold[TAU]);

                printf("\n");
            }

            printf("E_int (total)      = %e \n", E_int_from_cons(cp->cons));
            if(cp->multiphase)
            {
                printf("E_int (hot)        = %e \n", E_int_from_cons(cp->cons_hot));
                printf("E_int (cold)       = %e \n", E_int_from_cons(cp->cons_cold));
            }

            printf("\n");

            printf("v (total)          = %e \n", cp->cons[SRR]      / cp->cons[DDD]);
            if(cp->multiphase)
            {
                printf("v (hot)            = %e \n", cp->cons_hot[SRR]  / cp->cons_hot[DDD]);
                printf("v (cold)           = %e \n", cp->cons_cold[SRR] / cp->cons_cold[DDD]);
            }

            printf("\n");
        }


        double E_int_before_hot_L;  
        double E_int_before_cold_L; 
        double E_int_before_L;    

        double E_kin_before_hot_L;  
        double E_kin_before_cold_L; 
        double E_kin_before_L;  

        double M_before_hot_L;  
        double M_before_cold_L; 
        double M_before_L;    

        double Z_before_hot_L;  
        double Z_before_cold_L; 
        double Z_before_L;   


        if (c->multiphase)
        {
            E_int_before_hot_L  = E_int_from_cons(c->cons_hot);
            E_int_before_cold_L = E_int_from_cons(c->cons_cold);
            E_int_before_L      = E_int_from_cons(c->cons);

            E_kin_before_hot_L  = c->cons_hot[TAU]  - E_int_before_hot_L;
            E_kin_before_cold_L = c->cons_cold[TAU] - E_int_before_cold_L;
            E_kin_before_L      = c->cons[TAU]      - E_int_before_L;

            M_before_hot_L  = c->cons_hot[DDD];
            M_before_cold_L = c->cons_cold[DDD];
            M_before_L      = c->cons[DDD];

            Z_before_hot_L  = c->cons_hot[ZZZ];
            Z_before_cold_L = c->cons_cold[ZZZ];
            Z_before_L      = c->cons[ZZZ];
        }


        //Remove Zone at iS+1
        c->dr   += cp->dr;
        c->riph  = cp->riph;

        // if(c->multiphase || cp->multiphase)
        // {
        //     if(c->multiphase == 0)
        //     {
        //         printf("Left cell *not* multiphase\n");
        //         for( int q=0 ; q<NUM_Q ; ++q )
        //         {
        //             c->cons_cold[q]   = c->cons[q];
        //             c->RKcons_cold[q] = c->RKcons[q];
        //         }
        //     }
        //     else{
        //         printf("Left cell multiphase\n");
        //     }

        //     if(cp->multiphase)
        //     {
        //         printf("Right cell multiphase\n");
        //         for( int q=0 ; q<NUM_Q ; ++q )
        //         {
        //             c->cons_cold[q]   += cp->cons_cold[q];
        //             c->RKcons_cold[q] += cp->RKcons_cold[q];

        //             c->cons_hot[q]    += cp->cons_hot[q];
        //             c->RKcons_hot[q]  += cp->RKcons_hot[q];
        //         }
        //     }
        //     else
        //     {
        //         printf("Right cell *not* multiphase\n");

        //         for( int q=0 ; q<NUM_Q ; ++q )
        //         {
        //             c->cons_cold[q]   += cp->cons[q];
        //             c->RKcons_cold[q] += cp->RKcons[q];
        //         }
        //     }

        //     c->multiphase = 1;
        // }

        for( int q=0 ; q<NUM_Q ; ++q )
        {
            c->cons[q]   += cp->cons[q];
            c->RKcons[q] += cp->RKcons[q];
        }



        const double gamma = theDomain->theParList.Adiabatic_Index;
        const double rp = c->riph;
        const double rm = rp - c->dr;
        const double dV = get_dV( rp , rm );
        cons2prim( c->cons , c->prim , dV );
        if (c->multiphase)
        {
            const double x_hot_before  = M_before_hot_L  / M_before_L;
            const double x_cold_before = M_before_cold_L / M_before_L;

            const double y_hot_before  = E_int_before_hot_L  / E_int_before_L;
            const double y_cold_before = E_int_before_cold_L / E_int_before_L;

            const double z_hot_before  = Z_before_hot_L  / M_before_hot_L;
            const double z_cold_before = Z_before_cold_L / M_before_cold_L;



            const double E_int_after_L   = E_int_from_cons(c->cons);
            // const double E_kin_after_L   = c->cons[TAU] - E_int_after_L;
            const double E_kin_after_L   = E_kin_from_cons(c->cons);
            const double M_after_L       = c->cons[DDD];
            const double Z_after_L       = c->cons[ZZZ];

            const double dE_int_L   = E_int_after_L - E_int_before_L;
            const double dE_kin_L   = E_kin_after_L - E_kin_before_L;
            const double dM_L       = M_after_L - M_before_L;
            const double dZ_L       = Z_after_L - Z_before_L;

            printf(" --------- Transferring -------- \n");

            printf("dM_L      = %e \n", dM_L);
            printf("dZ_L      = %e \n", dZ_L);

            printf("dE_int_L      = %e \n", dE_int_L);
            printf("dE_kin_L      = %e \n", dE_kin_L);
            printf("\n");


            c->cons_hot[DDD]  += dM_L * x_hot_before;
            c->cons_cold[DDD] += dM_L * x_cold_before;

            c->cons_hot[SRR]  = c->cons[SRR] * (c->cons_hot[DDD]  / c->cons[DDD]);
            c->cons_cold[SRR] = c->cons[SRR] * (c->cons_cold[DDD] / c->cons[DDD]);


            c->cons_hot[TAU]  += (dE_int_L * y_hot_before)  + (dE_kin_L * x_hot_before);
            c->cons_cold[TAU] += (dE_int_L * y_cold_before) + (dE_kin_L * x_cold_before);

            c->cons_hot[ZZZ]  += dZ_L * x_hot_before * z_hot_before;
            c->cons_cold[ZZZ] += dZ_L * x_cold_before * z_cold_before;

            calc_multiphase_prim( c, c->prim_hot , c->prim_cold ,
                                     &(c->V_hot) , &(c->V_cold) );   


        }
        c->E_int_old  = c->prim[PPP] * dV / (gamma-1);
        c->dV_old     = dV;



        if(c->multiphase)
        {
            printf("------ After ------\n");

            printf("\n");

            printf("c->cons[DDD]       = %e \n", c->cons[DDD]);
            printf("c->cons[SRR]       = %e \n", c->cons[SRR]);
            printf("c->cons[TAU]       = %e \n", c->cons[TAU]);

            printf("\n");

            printf("c->cons_hot[DDD]   = %e \n", c->cons_hot[DDD]);
            printf("c->cons_hot[SRR]   = %e \n", c->cons_hot[SRR]);
            printf("c->cons_hot[TAU]   = %e \n", c->cons_hot[TAU]);

            printf("\n");

            printf("c->cons_cold[DDD]  = %e \n", c->cons_cold[DDD]);
            printf("c->cons_cold[SRR]  = %e \n", c->cons_cold[SRR]);
            printf("c->cons_cold[TAU]  = %e \n", c->cons_cold[TAU]);

            printf("\n");

            printf("E_int (total)      = %e \n", E_int_from_cons(c->cons));
            printf("E_int (hot)        = %e \n", E_int_from_cons(c->cons_hot));
            printf("E_int (cold)       = %e \n", E_int_from_cons(c->cons_cold));

            printf("\n");

            printf("v (total)          = %e \n", c->cons[SRR]      / c->cons[DDD]);
            printf("v (hot)            = %e \n", c->cons_hot[SRR]  / c->cons_hot[DDD]);
            printf("v (cold)           = %e \n", c->cons_cold[SRR] / c->cons_cold[DDD]);

            printf("\n");
        }


       // ----------- Post-conditions ------------- //
        if (c->multiphase)
        {
            verify_multiphase_conditions(c, "AMR()", "post");

            if (E_int_from_cons(c->cons_hot) < 0)
            {
                printf("------ERROR in AMR() postconditions------- \n");
                printf("hot gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_hot) = %e \n", E_int_from_cons(c->cons_hot));
                assert( E_int_from_cons(c->cons_hot) > 0 );
            }

            if (E_int_from_cons(c->cons_cold) < 0)
            {
                printf("------ERROR in AMR() postconditions------- \n");
                printf("cold gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_cold) = %e \n", E_int_from_cons(c->cons_cold));
                assert( E_int_from_cons(c->cons_cold) > 0 );
            }
        }



        //Shift Memory
        int blocksize = Nr-iSp-1;
        memmove( theCells+iSp , theCells+iSp+1 , blocksize*sizeof(struct cell) );
        theDomain->Nr -= 1;
        Nr = theDomain->Nr;
        theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
        theCells = theDomain->theCells;
        if( iS < iL ) iL--;

    }

    if( L>MaxLong )
    {
        if (iL == 1)
        {
            // printf("FORGE! iL = %d\n",iL);         
        }
        printf("FORGE! iL = %d\n",iL);

        theDomain->Nr += 1;
        Nr = theDomain->Nr;
        theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
        theCells = theDomain->theCells;
        int blocksize = Nr-iL-1;
        memmove( theCells+iL+1 , theCells+iL , blocksize*sizeof(struct cell) );

        struct cell * c  = &(theCells[iL]);
        struct cell * cp = &(theCells[iL+1]);

        double rp = c->riph;
        double rm = rp - c->dr;
        double r0 = pow( .5*(rp*rp*rp+rm*rm*rm) , 1./3. );

        double dV_orig = get_dV(rp, rm);
        double prim_tmp[NUM_Q];
        cons2prim ( c->cons , prim_tmp , dV_orig, true);

        c->riph  = r0;
        c->dr    = r0-rm;
        cp->dr   = rp-r0; // cp->riph already set at rp

        for( int q=0 ; q<NUM_Q ; ++q )
        {
            c->cons[q]    *= .5;
            c->RKcons[q]  *= .5;
            cp->cons[q]   *= .5;
            cp->RKcons[q] *= .5;

            c->cons_hot[q]    *= .5;
            c->RKcons_hot[q]  *= .5;
            cp->cons_hot[q]   *= .5;
            cp->RKcons_hot[q] *= .5;

            c->cons_cold[q]    *= .5;
            c->RKcons_cold[q]  *= .5;
            cp->cons_cold[q]   *= .5;
            cp->RKcons_cold[q] *= .5;

        }

        const double gamma = theDomain->theParList.Adiabatic_Index;
        double dV = get_dV( r0 , rm );
        cons2prim( c->cons , c->prim , dV , true);
        if (c->multiphase)
        {
            calc_multiphase_prim( c, c->prim_hot , c->prim_cold ,
                                     &(c->V_hot) , &(c->V_cold) );   
        }
        c->E_int_old  = cp->prim[PPP] * dV / (gamma-1);
        c->dV_old     = dV;

        dV = get_dV( rp , r0 );
        cons2prim( cp->cons , cp->prim , dV , true);
        if (cp->multiphase)
        {
            calc_multiphase_prim( cp, cp->prim_hot , cp->prim_cold ,
                                      &(cp->V_hot) , &(cp->V_cold) );   
        }
        cp->E_int_old  = cp->prim[PPP] * dV / (gamma-1);
        cp->dV_old     = dV;






        // ======== Verify post-conditions ========= //
        #ifndef NDEBUG
        // equal primitives + conservatives between two split cells
        for ( int q=0 ; q<NUM_Q ; ++q)
        {
            double tol = 1e-2;
            if ( std::abs((dV_orig/2 - dV) / dV) > tol)
            {
                printf("-----ERROR in AMR (forge) ------- \n");
                printf("just split cell: %d \n", iL);
                printf("expected dV_orig/2 == dV \n");
                printf("instead found: \n");
                printf("dV_orig/2  = %e \n", dV_orig/2);
                printf("dV         = %e \n", dV);
                printf("fractional error : %e \n", (dV_orig/2 - dV) / dV);
                assert(0);
            }
            if ( std::abs((c->prim[q] - prim_tmp[q]) / prim_tmp[q]) > tol)
            {
                printf("-----ERROR in AMR (forge) ------- \n");
                printf("just split cell: %d \n", iL);
                printf("expected c->prim[%d] == prim_tmp[%d] \n", q, q);
                printf("instead found: \n");
                printf("c->prim[%d]  = %e \n", q, c->prim[q]);
                printf("prim_tmp[%d] = %e \n", q, prim_tmp[q]);
                printf("fractional error : %e \n", (c->prim[q] - prim_tmp[q]) / prim_tmp[q]);
                assert(0);
            }
            if ( std::abs((cp->prim[q] - prim_tmp[q]) / prim_tmp[q]) > tol)
            {
                printf("-----ERROR in AMR (forge) ------- \n");
                printf("just split cell: %d \n", iL);
                printf("expected cp->prim[%d] == prim_tmp[%d] \n", q, q);
                printf("instead found: \n");
                printf("cp->prim[%d] = %e \n", q, cp->prim[q]);
                printf("prim_tmp[%d] = %e \n", q, prim_tmp[q]);
                printf("fractional error : %e \n", (cp->prim[q] - prim_tmp[q]) / prim_tmp[q]);
                assert(0);
            }
            if ( std::abs((cp->prim[q] - c->prim[q]) / c->prim[q]) > tol)
            {
                printf("-----ERROR in AMR (forge) ------- \n");
                printf("just split cell: %d \n", iL);
                printf("expected cp->prim[%d] == c->prim[%d] \n", q, q);
                printf("instead found: \n");
                printf("cp->prim[%d] = %e \n", q, cp->prim[q]);
                printf("c ->prim[%d] = %e \n", q, c ->prim[q]);
                printf("fractional error : %e \n", (cp->prim[q] - c->prim[q]) / c->prim[q]);
                assert(0);
            }
            if ( std::abs((cp->cons[q] - c->cons[q]) / c->cons[q]) > tol)
            {
                printf("-----ERROR in AMR (forge) ------- \n");
                printf("just split cell: %d \n", iL);
                printf("expected cp->cons[%d] == c->cons[%d] \n", q, q);
                printf("instead found: \n");
                printf("cp->cons[%d] = %e \n", q, cp->cons[q]);
                printf("c ->cons[%d] = %e \n", q, c ->cons[q]);
                printf("fractional error : %e \n", (cp->cons[q] - c->cons[q]) / c->cons[q]);
                assert(0);
            }
        }
        #endif
    }
}

void check_multiphase(struct domain * theDomain)
{
    // ============================================= //
    //
    //  For each cell: if it's multiphase check, check if it should be
    //  reverted back to single phase. If so, revert it.
    //  
    //  For multiphase cells that remain multiphase, enforce equal velocities
    //  in the cold + hot phases, based off the overall (single phase) velocity.
    //  In theory this should already be true, but numeric errors might have
    //  accumulated.
    //  
    //
    // ============================================= // 

    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;


    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        struct cell * c = &(theCells[i]);

        if(c->multiphase)
        {

            if (E_int_from_cons(c->cons_hot) < 0)
            {
                printf("------ERROR in check_multiphase() pre-conditions------- \n");
                printf("hot gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_hot) = %e \n", E_int_from_cons(c->cons_hot));
                assert( E_int_from_cons(c->cons_hot) > 0 );
            }

            if (E_int_from_cons(c->cons_cold) < 0)
            {
                printf("------ERROR in check_multiphase() pre-conditions------- \n");
                printf("cold gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_cold) = %e \n", E_int_from_cons(c->cons_cold));
                assert( E_int_from_cons(c->cons_cold) > 0 );
            }

            printf("i = %d \n", i);
            printf("--------------- check_multiphase - before SRR fix ---------- \n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]); 
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]); 

            printf("c->cons_hot[SRR]  = %e \n", c->cons_hot[SRR]);
            printf("c->cons_cold[SRR] = %e \n", c->cons_cold[SRR]); 
            printf("c->cons[SRR]      = %e \n", c->cons[SRR]); 


            // ensure that velocities are matched
            c->cons_hot[SRR]  = c->cons[SRR] * c->cons_hot[DDD]  / c->cons[DDD];
            c->cons_cold[SRR] = c->cons[SRR] * c->cons_cold[DDD] / c->cons[DDD];

            printf("--------------- check_multiphase - after SRR fix ---------- \n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]); 
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]); 

            printf("c->cons_hot[SRR]  = %e \n", c->cons_hot[SRR]);
            printf("c->cons_cold[SRR] = %e \n", c->cons_cold[SRR]); 
            printf("c->cons[SRR]      = %e \n", c->cons[SRR]); 

            verify_multiphase_conditions(c, "check_multiphase()", "mid");

            if (E_int_from_cons(c->cons_hot) < 0)
            {
                printf("------ERROR in check_multiphase() midconditions------- \n");
                printf("hot gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_hot) = %e \n", E_int_from_cons(c->cons_hot));
                assert( E_int_from_cons(c->cons_hot) > 0 );
            }

            if (E_int_from_cons(c->cons_cold) < 0)
            {
                printf("------ERROR in check_multiphase() midconditions------- \n");
                printf("cold gas internal energy less than 0\n");
                printf("E_int_from_cons(c->cons_cold) = %e \n", E_int_from_cons(c->cons_cold));
                assert( E_int_from_cons(c->cons_cold) > 0 );
            }


            const double T_hot  = calc_T(c->prim_hot);
            const double T_cold = calc_T(c->prim_cold);

            if ( (T_hot < 1e5) || (T_cold > 1e5) || ( c->cons_hot[DDD]/c->cons_cold[DDD]) > 1e5 ) 
            {
                // revert to single phase
                c->multiphase=0;
                for ( int q=0 ; q<NUM_Q ; ++q)
                {
                    c->prim_hot[q]  = 0;
                    c->prim_cold[q] = 0;
                    c->cons_hot[q]  = 0;
                    c->cons_cold[q] = 0;
                    c->RKcons_hot[q]  = 0;
                    c->RKcons_cold[q] = 0;
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
        }
    }
}


