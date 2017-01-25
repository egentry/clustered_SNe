
#include <cmath>
#include <assert.h>
#include <stdlib.h>
extern "C" {
#include <grackle.h>
}

#include "../structure.H"
#include "../cooling.H"
#include "../misc.H" // get_mean_molecular_weight
#include "../geometry.H" // get_dV
#include "../constants.H" // m_proton, k_boltzmann
#include "euler.H"

static double GAMMA_LAW = 0.0;
static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;

static double H_0 = 0.0;
static double H_1 = 0.0;


void setHydroParams( const struct domain * theDomain )
{
    GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
    RHO_FLOOR = theDomain->theParList.Density_Floor;
    PRE_FLOOR = theDomain->theParList.Pressure_Floor;

    H_0 = theDomain->theParList.H_0;
    H_1 = theDomain->theParList.H_1;
}

void prim2cons( const double * prim , double * cons , const double dV )
{
    // prim2cons shouldn't overwrite the cons array of any cell,
    // besides when setting up initial conditions.

    const double rho  = prim[RHO];
    const double P    = prim[PPP];
    const double vr   = prim[VRR];
    const double v2   = vr*vr;
    const double gam  = GAMMA_LAW;
    const double rhoe = P/(gam-1.);

    cons[DDD] = rho*dV;
    cons[SRR] = rho*vr*dV;
    cons[TAU] = (.5*rho*v2 + rhoe)*dV;

    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        cons[q] = cons[DDD]*prim[q];
    }
}

void cons2prim( const double * cons , double * prim , 
                const double dV , const bool verbose)
{
    // prim2cons shouldn't overwrite the cons array of any cell,
    // besides when setting up initial conditions.
    // Any of the RHO_FLOOR, PRE_FLOOR adjustments below
    // shouldn't actually change the cons values 

    // // E    :    total energy / unit volume
    // // e    : internal energy / unit mass
    // // rhoe : internal energy / unit volume
    // // P    : pressure
    double rho = cons[DDD]/dV;
    const double Sr  = cons[SRR]/dV;
    const double E   = cons[TAU]/dV;

    const double vr   = Sr/rho;
    const double v2   = vr*vr;
    const double rhoe = E - .5*rho*v2;
    const double gam  = GAMMA_LAW;
    double P    = (gam-1.)*rhoe;

    if( rho<RHO_FLOOR ) 
    {
        // any changes here just affect flux calculations,
        // NOT the actual mass

        printf("------ ERROR in cons2prim() -- RHO_FLOOR ------- \n");
        printf("rho = %e \n", rho );
        printf("Expected rho > %e \n", RHO_FLOOR);
        printf("dV = %e \n", dV);
        rho=RHO_FLOOR; 

        // assert(rho>RHO_FLOOR);
    }
    if( P < PRE_FLOOR )
    {
        // any changes here just affect flux calculations,
        // NOT the actual internal energy
        if (verbose)
        {
            printf("------ ERROR in cons2prim()------- \n");
            printf("pressure should be above pressure floor! \n");
            printf("pressure       = %e \n", P);
            printf("pressure floor = %e \n", PRE_FLOOR);
            printf("dV  = %e \n", dV);
            printf("rho = %e \n", rho);
            printf("vr  = %e \n", vr);
        }
        // assert(P  > 0);

        P = PRE_FLOOR;
    } 

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[VRR] = vr;

    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        prim[q] = cons[q]/cons[DDD];
    }

}

void calc_multiphase_prim(const struct cell * c, 
    double * prim_hot, double * prim_cold,
    double * V_hot,    double * V_cold)
{
    // input: cell c
    // outputs: prim_hot, prim_cold array, V_hot, V_cold pointer to double (non-array)
    //
    // assumes c->prim, c->cons_hot, c->cons_cold already up to date
    //
    // has no side effects, besides prim_hot, prim_cold, V_hot, V_cold input pointers
    if( c->cons_hot[DDD]<=0 ) 
    {
        printf("------ ERROR in calc_multiphase_prim() preconditions ------- \n");
        printf("Hot mass less than 0; \n");

        printf("c->cons_hot[DDD]                     = %e \n", c->cons_hot[DDD]);
        printf("c->cons_cold[DDD]                    = %e \n", c->cons_cold[DDD]);
        printf("c->cons_hot[DDD] + c->cons_cold[DDD] = %e \n", c->cons_hot[DDD]+c->cons_cold[DDD]);
        printf("c->cons[DDD]                         = %e \n", c->cons[DDD]);

        assert( c->cons_hot[DDD]>0 );
    }

    if( c->cons_cold[DDD]<=0 ) 
    {
        printf("------ ERROR in calc_multiphase_prim() preconditions ------- \n");
        printf("Cold mass less than 0; \n");

        printf("c->cons_hot[DDD]                     = %e \n", c->cons_hot[DDD]);
        printf("c->cons_cold[DDD]                    = %e \n", c->cons_cold[DDD]);
        printf("c->cons_hot[DDD] + c->cons_cold[DDD] = %e \n", c->cons_hot[DDD]+c->cons_cold[DDD]);
        printf("c->cons[DDD]                         = %e \n", c->cons[DDD]);

        assert( c->cons_cold[DDD]>0 );
    }

    if( E_int_from_cons(c->cons_hot) <=0 ) 
    {
        printf("------ ERROR in calc_multiphase_prim() preconditions ------- \n");
        printf("Hot internal energy less than 0; \n");
        printf("c->cons_hot[TAU]               = %e \n", c->cons_hot[TAU]);
        printf("E_int_from_cons(c->cons_hot)   = %e \n", E_int_from_cons(c->cons_hot));

        printf("c->cons[TAU]                   = %e \n", c->cons[TAU]);
        printf("E_int_from_cons(c->cons)       = %e \n", E_int_from_cons(c->cons));

        assert( E_int_from_cons(c->cons_hot) >0 );
    }

    if( E_int_from_cons(c->cons_cold) <=0 ) 
    {
        printf("------ ERROR in calc_multiphase_prim() preconditions ------- \n");
        printf("Cold internal energy less than 0; \n");
        printf("c->cons_cold[TAU]              = %e \n", c->cons_cold[TAU]);
        printf("E_int_from_cons(c->cons_cold)  = %e \n", E_int_from_cons(c->cons_cold));

        printf("c->cons[TAU]                   = %e \n", c->cons[TAU]);
        printf("E_int_from_cons(c->cons)       = %e \n", E_int_from_cons(c->cons));

        assert( E_int_from_cons(c->cons_cold) >0 );
    }


    prim_hot[VRR]  = c->prim[VRR];
    prim_cold[VRR] = c->prim[VRR];

    prim_hot[ZZZ]  = c->cons_hot[ZZZ]  / c->cons_hot[DDD];
    prim_cold[ZZZ] = c->cons_cold[ZZZ] / c->cons_cold[DDD];

    prim_hot[PPP]  = c->prim[PPP];
    prim_cold[PPP] = c->prim[PPP];

    const double E_int_hot  = E_int_from_cons(c->cons_hot);
    const double E_int_cold = E_int_from_cons(c->cons_cold);

    prim_hot[RHO]  = prim_hot[PPP]  * c->cons_hot[DDD]  / (GAMMA_LAW-1) / E_int_hot;
    prim_cold[RHO] = prim_cold[PPP] * c->cons_cold[DDD] / (GAMMA_LAW-1) / E_int_cold;

    *V_hot  = c->cons_hot[DDD]  / prim_hot[RHO];
    *V_cold = c->cons_cold[DDD] / prim_cold[RHO];

    if( prim_hot[RHO]<RHO_FLOOR ) 
    {

        printf("------ ERROR in calc_multiphase_prim() -- RHO_FLOOR (post-conditions) ------- \n");
        printf("prim_hot[RHO]    = %e \n", prim_hot[RHO] );
        printf("prim_hot[PPP]    = %e \n", prim_hot[PPP] );
        printf("c->cons_hot[DDD] = %e \n", c->cons_hot[DDD]);

        printf("Expected rho > %e \n", RHO_FLOOR);

        assert(prim_hot[RHO]>RHO_FLOOR);
    }

    if( prim_cold[RHO]<RHO_FLOOR ) 
    {

        printf("------ ERROR in calc_multiphase_prim() -- RHO_FLOOR (post-conditions) ------- \n");
        printf("prim_cold[RHO]    = %e \n", prim_cold[RHO] );
        printf("prim_cold[PPP]    = %e \n", prim_cold[PPP] );
        printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
        printf("Expected rho > %e \n", RHO_FLOOR);

        assert(prim_cold[RHO]>RHO_FLOOR);
    }

}

double E_int_from_cons( const double * cons )
{

    const double mass  = cons[DDD];
    const double vr    = cons[SRR] / mass;
    const double E_int = cons[TAU] - .5*mass*vr*vr;

    return E_int;
}

double E_int_from_prim( const struct cell * c , 
                        const double gamma )
{

    const double rp = c->riph;
    const double rm = rp - c->dr;
    const double dV = get_dV( rp , rm );
    const double E_int = c->prim[PPP] * dV / (gamma-1);

    assert(E_int > 0);

    return E_int;
}

void getUstar( const double * prim ,  double * Ustar ,
               const double Sk , const double Ss )
{

    const double rho  = prim[RHO];
    const double vr   = prim[VRR];
    const double P    = prim[PPP];
    const double v2   = vr*vr;

    const double gam  = GAMMA_LAW;

    const double rhoe = P/(gam-1.);

    const double rhostar =  rho*(Sk - vr)/(Sk - Ss);
    const double Pstar   =    P*(Ss - vr)/(Sk - Ss);
    const double Us      = rhoe*(Sk - vr)/(Sk - Ss);

    Ustar[DDD] = rhostar;
    Ustar[SRR] = rhostar*Ss;
    Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss*(Ss - vr) + Pstar;

    // ======== Verify post-conditions ========= //
    // require finite fluxes
    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        Ustar[q] = prim[q]*Ustar[DDD];
    }

    #ifndef NDEBUG
    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        if( !std::isfinite(Ustar[q]) )
        {
            printf("Ustar[%d] = %e in 'getUstar()' \n", q, Ustar[q]);
            printf("prim[%d]  = %e in 'getUstar()' \n", q, prim[q]);
            printf("Sk        = %20.10le in 'getUstar()' \n", Sk);
            printf("Ss        = %20.10le in 'getUstar()' \n", Ss);
            assert( std::isfinite(Ustar[q]) );
        }
    }
    #endif

}

void flux( const double * prim , double * flux )
{
    const double rho  = prim[RHO];
    const double P    = prim[PPP];
    const double vr   = prim[VRR];
    const double v2   = vr*vr;
    const double gam  = GAMMA_LAW;
    const double rhoe = P/(gam-1.);

    flux[DDD] = rho*vr;
    flux[SRR] = rho*vr*vr + P;
    flux[TAU] = (.5*rho*v2 + rhoe + P)*vr;

    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        flux[q] = flux[DDD]*prim[q];
    }
}



void source( struct cell * c, 
             const double rp , const double rm ,
             const double dV , const double dt ,
             const double R_shock, Cooling * cooling )
{
    if (c->multiphase)
    {
        // printf("------- In source (start) ---------- \n");
        // printf("E_int (cold)      = %e \n", E_int_from_cons(c->cons_cold));
        // printf("E_int (hot)       = %e \n", E_int_from_cons(c->cons_hot));
        // printf("E_int (total)     = %e \n", E_int_from_cons(c->cons));

        // assert(E_int_from_cons(c->cons_cold) > 0);
        // assert(E_int_from_cons(c->cons_hot) > 0);
        assert(E_int_from_cons(c->cons) > 0);   
    }


    const double P = c->prim[PPP]; // pressure
    const double r = .5 * (rp + rm);

    double dP = 0; // change in momentum, not pressure
    // 1st order contribution
    dP += 4 * M_PI * P * (std::pow(rp, 2.) - std::pow(rm, 2.)) * dt;
    // 2nd order contribution
    dP += 8*M_PI*c->grad[PPP]
                *(     (std::pow(rp,3.) - std::pow(rm,3.))/3. 
                   + r*(std::pow(rp,2.) - std::pow(rm,2.))/2. );

    c->cons[SRR] += dP;
    if (c->multiphase)
    {
        // c->cons_hot[SRR]  += dP * (c->cons_hot[DDD]  / c->cons[DDD]);
        // c->cons_cold[SRR] += dP * (c->cons_cold[DDD] / c->cons[DDD]);

        // c->cons_hot[SRR]  += dP * c->x_hot;
        // c->cons_cold[SRR] += dP * c->x_cold;

        // c->cons_hot[SRR]  = c->cons[SRR] * c->x_hot;
        // c->cons_cold[SRR] = c->cons[SRR] * c->x_cold;

        c->cons_hot[SRR]  = c->cons[SRR] * (c->cons_hot[DDD] / c->cons[DDD]);
        c->cons_cold[SRR] = c->cons[SRR] * (c->cons_cold[DDD] / c->cons[DDD]);

        // assert(E_int_from_cons(c->cons_cold) > 0);
        // assert(E_int_from_cons(c->cons_hot) > 0);
        // assert(E_int_from_cons(c->cons) > 0);
    }       

    // now decompose the change in energy, and apply it to subgrid components
    if (c->multiphase)
    {
        const double E_int = E_int_from_cons(c->cons);
        const double E_kin = c->cons[TAU] - E_int;

        const double dE_int = E_int - c->E_int_initial;
        const double dE_kin = E_kin - c->E_kin_initial;

        printf("------- In source (*before* apply decomposed energy) ---------- \n");
        printf("E_int (hot)       = %e \n", E_int_from_cons(c->cons_hot));
        printf("E_int (cold)      = %e \n", E_int_from_cons(c->cons_cold));
        printf("E_int (total)     = %e \n", E_int_from_cons(c->cons));

        c->cons_hot[TAU]  += (c->y_hot  * dE_int) + (c->x_hot  * dE_kin);      
        c->cons_cold[TAU] += (c->y_cold * dE_int) + (c->x_cold * dE_kin);

        printf("------- In source (*after* apply decomposed energy) ---------- \n");
        printf("E_int (hot)       = %e \n", E_int_from_cons(c->cons_hot));
        printf("E_int (cold)      = %e \n", E_int_from_cons(c->cons_cold));
        printf("E_int (total)     = %e \n", E_int_from_cons(c->cons));

            // printf("\n");


            // printf("c->cons_hot[DDD]     = %e \n", c->cons_hot[DDD]);
            // printf("c->cons_hot[SRR]     = %e \n", c->cons_hot[SRR]);
            // printf("c->cons_hot[TAU]     = %e \n", c->cons_hot[TAU]);
            // printf("c->cons_hot[ZZZ]     = %e \n", c->cons_hot[ZZZ]);

            // printf("\n");

            // printf("c->cons_cold[DDD]    = %e \n", c->cons_cold[DDD]);
            // printf("c->cons_cold[SRR]    = %e \n", c->cons_cold[SRR]);
            // printf("c->cons_cold[TAU]    = %e \n", c->cons_cold[TAU]);
            // printf("c->cons_cold[ZZZ]    = %e \n", c->cons_cold[ZZZ]);

            // printf("\n");

            // printf("c->cons[DDD]         = %e \n", c->cons[DDD]);
            // printf("c->cons[SRR]         = %e \n", c->cons[SRR]);
            // printf("c->cons[TAU]         = %e \n", c->cons[TAU]);
            // printf("c->cons[ZZZ]         = %e \n", c->cons[ZZZ]);

            // printf("\n");

            // printf("VRR (hot)         = %e \n", c->cons_hot[SRR]  / c->cons_hot[DDD]);
            // printf("VRR (cold)        = %e \n", c->cons_cold[SRR] / c->cons_cold[DDD]);
            // printf("VRR (total)       = %e \n", c->cons[SRR]      / c->cons[DDD]);


        assert(E_int_from_cons(c->cons_hot)  > 0); 
        assert(E_int_from_cons(c->cons_cold) > 0);
        assert(E_int_from_cons(c->cons)      > 0); 

    }


    if( cooling->with_cooling == true )
    {
        if (c->multiphase)
        {
            bool cached_cooling = false;

            // calc_multiphase_prim(c, c->prim_hot, c->prim_cold,
            //                         &(c->V_hot), &(c->V_cold) );

            printf("----- In source (multiphase cooling) ------- \n");

            printf("c->cons_hot[DDD]    = %e \n", c->cons_hot[DDD]);
            printf("c->cons_hot[SRR]    = %e \n", c->cons_hot[SRR]);
            printf("c->cons_hot[TAU]    = %e \n", c->cons_hot[TAU]);
            printf("c->cons_hot[ZZZ]    = %e \n", c->cons_hot[ZZZ]);

            printf("\n");

            printf("c->prim_hot[RHO]    = %e \n", c->prim_hot[RHO]);
            printf("c->prim_hot[VRR]    = %e \n", c->prim_hot[VRR]);
            printf("c->prim_hot[PPP]    = %e \n", c->prim_hot[PPP]);
            printf("c->prim_hot[ZZZ]    = %e \n", c->prim_hot[ZZZ]);

            printf("\n");

            printf("c->cons_cold[DDD]    = %e \n", c->cons_cold[DDD]);
            printf("c->cons_cold[SRR]    = %e \n", c->cons_cold[SRR]);
            printf("c->cons_cold[TAU]    = %e \n", c->cons_cold[TAU]);
            printf("c->cons_cold[ZZZ]    = %e \n", c->cons_cold[ZZZ]);

            printf("\n");

            printf("c->prim_cold[RHO]    = %e \n", c->prim_cold[RHO]);
            printf("c->prim_cold[VRR]    = %e \n", c->prim_cold[VRR]);
            printf("c->prim_cold[PPP]    = %e \n", c->prim_cold[PPP]);
            printf("c->prim_cold[ZZZ]    = %e \n", c->prim_cold[ZZZ]);

            printf("\n");

            printf("E_int (hot)          = %e \n", E_int_from_cons(c->cons_hot));
            printf("E_int (cold)         = %e \n", E_int_from_cons(c->cons_cold));
            printf("E_int (total)        = %e \n", E_int_from_cons(c->cons));

            printf("\n");

            calc_multiphase_prim( c, c->prim_hot , c->prim_cold ,
                                     &(c->V_hot) , &(c->V_cold) );    

            printf("----- In source (after recalculating multiphase_prim) ------- \n");


            printf("c->prim_hot[RHO]    = %e \n", c->prim_hot[RHO]);
            printf("c->prim_hot[VRR]    = %e \n", c->prim_hot[VRR]);
            printf("c->prim_hot[PPP]    = %e \n", c->prim_hot[PPP]);
            printf("c->prim_hot[ZZZ]    = %e \n", c->prim_hot[ZZZ]);

            printf("\n");

            printf("c->prim_cold[RHO]    = %e \n", c->prim_cold[RHO]);
            printf("c->prim_cold[VRR]    = %e \n", c->prim_cold[VRR]);
            printf("c->prim_cold[PPP]    = %e \n", c->prim_cold[PPP]);
            printf("c->prim_cold[ZZZ]    = %e \n", c->prim_cold[ZZZ]);

            printf("\n");           

            const double dE_cool_hot  = cooling->calc_cooling(c->prim_hot,  c->cons_hot,  dt , cached_cooling ) * c->V_hot;  
            const double dE_cool_cold = cooling->calc_cooling(c->prim_cold, c->cons_cold, dt , cached_cooling ) * c->V_cold;  

            printf("----- In source (after calculating multiphase cooling) ------- \n");


            printf("c->cons_hot[TAU]     = %e \n", c->cons_hot[TAU]);
            printf("c->cons_cold[TAU]    = %e \n", c->cons_cold[TAU]);
            printf("c->cons[TAU]         = %e \n", c->cons[TAU]);


            printf("\n");

            printf("dE_cool_hot          = %e \n", dE_cool_hot);
            printf("dE_cool_cold         = %e \n", dE_cool_cold);
            printf("dE_cool              = %e \n", dE_cool_hot + dE_cool_cold);

            printf("\n");

            printf("c->V_hot             = %e \n", c->V_hot);
            printf("c->V_cold            = %e \n", c->V_cold);
            printf("c->V_hot + c->V_cold = %e \n", c->V_hot + c->V_cold);
            printf("dV                   = %e \n", dV);

            c->dE_cool = dE_cool_hot + dE_cool_cold;
            c->cons_hot[TAU]  += dE_cool_hot;
            c->cons_cold[TAU] += dE_cool_cold;
            c->cons[TAU]      += c->dE_cool;


            printf("----- Exiting source (multiphase cooling) ------- \n");
            printf("E_int (hot)       = %e \n", E_int_from_cons(c->cons_hot));
            printf("E_int (cold)      = %e \n", E_int_from_cons(c->cons_cold));
            printf("E_int (total)     = %e \n", E_int_from_cons(c->cons));
            printf("-----\n\n");

            assert(E_int_from_cons(c->cons_cold) > 0);
            assert(E_int_from_cons(c->cons_hot) > 0);
            assert(E_int_from_cons(c->cons) > 0);

        }
        else
        {
            // never use cached cooling, so that I don't have to worry about
            // the previous cooling call being for a subgrid zone
            bool cached_cooling = false;
            // double dr = rp - rm;
            // if( r > (R_shock+(10*dr)) ) cached_cooling = true;

            c->dE_cool   = cooling->calc_cooling(c->prim, c->cons, dt , cached_cooling ) * dV;  
            c->cons[TAU] += c->dE_cool;   
        }

    }
    else
    {
        c->dE_cool = 0;
    }
}


double calc_T(const double * prim)
{

    // ============================================= //
    //
    //  Calculate temperature, assuming ideal gas of fixed mean molec. weight
    //
    //  Inputs:
    //    - prim      - array of primitive variables
    //
    //  Returns:
    //    - T         - Temperature [K]
    //
    //  Notes:
    //    - Assumes fixed mean molec. weight, set by cell metallicity
    //
    // ============================================= // 

    const double pressure = prim[PPP];
    const double density = prim[RHO];
    const double Z = prim[ZZZ];
    const double mu = get_mean_molecular_weight(Z);

    const double T = pressure / density * (mu * m_proton / k_boltzmann);

    return T;

}

void conduction( struct cell * cL , struct cell * cR, 
                 const double dA , const double dt )
{

    // ============================================= //
    //
    //  Calculate and add heat fluxes from *physical* conduction,
    //  using the equations of Keller et al. (2014), section 2 (eq. 4 & 6)
    //
    //  Inputs:
    //    - cL        - cell to the Left
    //    - cR        - cell to the Right
    //    - dA        - interface area between cells [cm^2]
    //    - dt        - timestep [s]
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - overwrites cons[DDD] (mass) and cons[TAU] (energy) 
    //      of both cells cL, cR (cL, cR on left and right side of interface, respectively)
    //    - if either cell is flagged as multiphase, it'll also overwrite
    //      cons_hot[DDD], cons_cold[DDD], cons_hot[TAU], cons_cold[TAU]
    //      of that cell
    //
    //  Notes:
    //    - Sign convention for dM_dt: always positive
    //
    // ============================================= // 


    // ======== Verify pre-conditions ========= //
    const double rel_tol = 1e-5; // relative tolerance for float comparison
    if (cL->multiphase)
    {
        if (cL->cons_hot[DDD] < 0)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[left] hot gas mass less than 0\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( cL->cons_hot[DDD] > 0 );
        }

        if (cL->cons_cold[DDD] < 0)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[left] cold gas mass less than 0\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( cL->cons_cold[DDD] > 0 );
        }

        if ( (cL->cons_hot[DDD] / cL->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[left] hot gas mass greater than total mass\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( (cL->cons_hot[DDD] / cL->cons[DDD]) < (1+rel_tol) );
        }

        if ( (cL->cons_cold[DDD] / cL->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[left] cold gas mass greater than total mass\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( (cL->cons_cold[DDD] / cL->cons[DDD]) < (1+rel_tol) );
        }

        if ( std::abs(1-( (cL->cons_cold[DDD] + cL->cons_hot[DDD])/cL->cons[DDD])) > rel_tol)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[left] cold mass + hot mass =/= total mass\n");
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert(  std::abs(1-( (cL->cons_cold[DDD] + cL->cons_hot[DDD])/cL->cons[DDD])) <= rel_tol);
        }




        if (cL->cons_hot[TAU] < 0)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[left] hot gas energy less than 0\n");
            printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

            assert( cL->cons_hot[TAU] > 0 );
        }

        if (cL->cons_cold[TAU] < 0)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[left] cold gas energy less than 0\n");
            printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

            assert( cL->cons_cold[TAU] > 0 );
        }

        // if (E_int_from_cons(cL->cons_hot) < 0)
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[left] hot gas *internal* energy less than 0\n");
        //     printf("[left] E_int (hot)  = %e \n", E_int_from_cons(cL->cons_hot));
        //     printf("cL->cons_hot[TAU]   = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons[TAU]       = %e \n", cL->cons[TAU]);

        //     assert( E_int_from_cons(cL->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(cL->cons_cold) < 0)
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[left] cold gas *internal* energy less than 0\n");
        //     printf("[left] E_int (cold)  = %e \n", E_int_from_cons(cL->cons_cold));
        //     printf("cL->cons_cold[TAU]   = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons[TAU]        = %e \n", cL->cons[TAU]);

        //     assert( E_int_from_cons(cL->cons_cold) > 0 );
        // }

        // if ( (cL->cons_hot[TAU] / cL->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[left] hot gas energy greater than total energy\n");
        //     printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

        //     assert( (cL->cons_hot[TAU] / cL->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( (cL->cons_cold[TAU] / cL->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[left] cold gas energy greater than total energy\n");
        //     printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

        //     assert( (cL->cons_cold[TAU] / cL->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( std::abs(1-( (cL->cons_cold[TAU] + cL->cons_hot[TAU])/cL->cons[TAU])) > rel_tol)
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[left] cold energy + hot energy =/= total energy\n");
        //     printf("cL->cons_cold[TAU]                      = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons_hot[TAU]                       = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons_cold[TAU] + cL->cons_hot[TAU]  = %e \n", cL->cons_cold[TAU] + cL->cons_hot[TAU]);
        //     printf("cL->cons[TAU]                           = %e \n", cL->cons[TAU]);
        //     printf("relative error  = %e \n", 1 - ( (cL->cons_cold[TAU] + cL->cons_hot[TAU]) / cL->cons[TAU]));

        //     assert(  std::abs(1-( (cL->cons_cold[TAU] + cL->cons_hot[TAU])/cL->cons[TAU])) <= rel_tol);
        // }
    }

    if (cR->multiphase)
    {
        if (cR->cons_hot[DDD] < 0)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[right] hot gas mass less than 0\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( cR->cons_hot[DDD] > 0 );
        }

        if (cR->cons_cold[DDD] < 0)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[right] cold gas mass less than 0\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( cR->cons_cold[DDD] > 0 );
        }

        if ( (cR->cons_hot[DDD] / cR->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[right] hot gas mass greater than total mass\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( (cR->cons_hot[DDD] / cR->cons[DDD]) < (1+rel_tol) );
        }

        if ( (cR->cons_cold[DDD] / cR->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[right] cold gas mass greater than total mass\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( (cR->cons_cold[DDD] / cR->cons[DDD]) < (1+rel_tol) );
        }

        if ( std::abs(1-( (cR->cons_cold[DDD] + cR->cons_hot[DDD])/cR->cons[DDD])) > rel_tol)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[right] cold mass + hot mass =/= total mass\n");
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert(  std::abs(1-( (cR->cons_cold[DDD] + cR->cons_hot[DDD])/cR->cons[DDD])) <= rel_tol);
        }




        if (cR->cons_hot[TAU] < 0)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[right] hot gas energy less than 0\n");
            printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
            printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
            printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

            assert( cR->cons_hot[TAU] > 0 );
        }

        if (cR->cons_cold[TAU] < 0)
        {
            printf("------ERROR in conduction() preconditions------- \n");
            printf("[right] cold gas energy less than 0\n");
            printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
            printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
            printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

            assert( cR->cons_cold[TAU] > 0 );
        }

        // if (E_int_from_cons(cR->cons_hot) < 0)
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[right] hot gas *internal* energy less than 0\n");
        //     printf("[right] E_int (hot) = %e \n", E_int_from_cons(cR->cons_hot));
        //     printf("cR->cons_hot[TAU]   = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons[TAU]       = %e \n", cR->cons[TAU]);

        //     assert( E_int_from_cons(cR->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(cR->cons_cold) < 0)
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[right] cold gas *internal* energy less than 0\n");
        //     printf("[right] E_int (cold) = %e \n", E_int_from_cons(cR->cons_cold));
        //     printf("cR->cons_cold[TAU]   = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons[TAU]        = %e \n", cR->cons[TAU]);

        //     assert( E_int_from_cons(cR->cons_cold) > 0 );
        // }

        // if ( (cR->cons_hot[TAU] / cR->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[right] hot gas energy greater than total energy\n");
        //     printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

        //     assert( (cR->cons_hot[TAU] / cR->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( (cR->cons_cold[TAU] / cR->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[right] cold gas energy greater than total energy\n");
        //     printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

        //     assert( (cR->cons_cold[TAU] / cR->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( std::abs(1-( (cR->cons_cold[TAU] + cR->cons_hot[TAU])/cR->cons[TAU])) > rel_tol)
        // {
        //     printf("------ERROR in conduction() preconditions------- \n");
        //     printf("[right] cold energy + hot energy =/= total energy\n");
        //     printf("cR->cons_cold[TAU]                      = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons_hot[TAU]                       = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons_cold[TAU] + cR->cons_hot[TAU]  = %e \n", cR->cons_cold[TAU] + cR->cons_hot[TAU]);
        //     printf("cR->cons[TAU]                           = %e \n", cR->cons[TAU]);
        //     printf("relative error  = %e \n", 1 - ( (cR->cons_cold[TAU] + cR->cons_hot[TAU]) / cR->cons[TAU]));

        //     assert(  std::abs(1-( (cR->cons_cold[TAU] + cR->cons_hot[TAU])/cR->cons[TAU])) <= rel_tol);
        // }
    }


    double prim_L[NUM_Q];
    double prim_R[NUM_Q];

    const double dr_L = .5*cL->dr;
    const double dr_R = .5*cR->dr;

    const double dx = dr_L + dr_R;

    const double T_L = calc_T(cL->prim);
    const double T_R = calc_T(cR->prim);

    if (std::abs(1 - (T_L/T_R)) < 1e4)
    {
        // effectively no temperature gradient; don't apply conduction
        return;
    }


    for( int q=0 ; q<NUM_Q ; ++q )
    {
        // extrapolate to boundary between cells
        prim_L[q] = cL->prim[q] + cL->grad[q]*dr_L;
        prim_R[q] = cR->prim[q] - cR->grad[q]*dr_R;
    }

    // calculate *UNsaturated* flux

    const double mu_L = get_mean_molecular_weight(prim_L[ZZZ]);
    const double mu_R = get_mean_molecular_weight(prim_R[ZZZ]);
    const double mu = (mu_L + mu_R) / 2.;

    const double dM_dt_unsat = (4.*M_PI * mu * m_proton) / (25. * k_boltzmann) * kappa_0
        * std::abs(std::pow(T_R, 5./2.) - std::pow(T_L, 5./2.)) / dx * dA; 

    // calculate *saturated* flux
    const double c_s_L = std::sqrt(std::abs(GAMMA_LAW*prim_L[PPP]/prim_L[RHO]));
    const double c_s_R = std::sqrt(std::abs(GAMMA_LAW*prim_R[PPP]/prim_R[RHO]));


    const double phi_s = 1.1; // see note near Eq. 8 of Cowie & McKee 1977
    double dM_dt_sat;
    if (T_L > T_R)
    {
        // left cell hotter -- use the [colder] right cell's sound speed + pressure
        const double q = 5 * phi_s * c_s_R * prim_R[PPP];
        dM_dt_sat = q * dA / (c_s_R*c_s_R); // it's the colder cell's energy which is transported
    }
    else
    {
        // right cell hotter -- use the [colder] left cell's sound speed + pressure
        const double q = 5 * phi_s * c_s_L * prim_L[PPP];
        dM_dt_sat = q * dA / (c_s_L*c_s_L);
    }

    const double dM_dt = std::min(dM_dt_unsat, dM_dt_sat);
    const double dM = dM_dt * dt;


    // to do: also transfer the correct amount of kinetic energy!!!


    double x; // mass fraction of transfered flux, relative to initial cell
    double dP; // transfered momentum
    double dE; // transfered energy
    double dM_Z; // transfered mass of metals
    int sgn; // sign; 1 if mass is transfered right, -1 if mass transfered left
    if (T_L > T_R)
    {
        sgn = -1;
        x = dM / cR->cons[DDD];
        dP = x * cR->cons[SRR];
        dE = (dM * c_s_R * c_s_R) + (dP*dP / (2*dM)); // internal + kinetic -- should I just transfer x * cons[TAU]?
        dM_Z = x * cR->cons[ZZZ];
    }
    else
    {
        sgn = 1;
        x = dM / cL->cons[DDD];
        dP = x * cL->cons[SRR];
        dE = (dM * c_s_L * c_s_L) + (dP*dP / (2*dM)); // internal + kinetic -- should I just transfer x * cons[TAU]?
        dM_Z = x * cL->cons[ZZZ];
    }


    // before applying fluxes, get multiphase mass + energy fractions
    // NOTE: in theory x_hot_L + x_cold_L = 1, but this gets tricky when they are
    // vastly different magnitudes. Rather than assuming x_cold_L = 1 - x_hot_L,
    // we'll calculate both x_hot_L and x_cold_L seperately, just in case.

    // double x_hot_L  = 0; // mass fraction of hot gas in left  cell, if multiphase
    // double x_hot_R  = 0; // mass fraction of hot gas in right cell, if multiphase

    // double y_hot_L  = 0; // energy fraction of hot gas in left  cell, if multiphase
    // double y_hot_R  = 0; // energy fraction of hot gas in right cell, if multiphase

    // double x_cold_L = 0; // mass fraction of cold gas in left  cell, if multiphase
    // double x_cold_R = 0; // mass fraction of cold gas in right cell, if multiphase

    // double y_cold_L = 0; // energy fraction of cold gas in left  cell, if multiphase
    // double y_cold_R = 0; // energy fraction of cold gas in right cell, if multiphase



    // if (cL->multiphase)
    // {
    //     x_hot_L  = cL->cons_hot[DDD] / cL->cons[DDD];
    //     y_hot_L  = E_int_from_cons(cL->cons_hot) / E_int_from_cons(cL->cons);

    //     x_cold_L = cL->cons_cold[DDD] / cL->cons[DDD];
    //     y_cold_L = E_int_from_cons(cL->cons_cold) / E_int_from_cons(cL->cons);

    //         if( (x_hot_L > 1 + rel_tol)  || (x_hot_L < 0))
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("x_hot_L not within (0,1 + rel_tol) \n");
    //             printf("x_hot_L = %e \n", x_hot_L);
    //             printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
    //             printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
    //             printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

    //             assert( (x_hot_L > 0) && (x_hot_L < 1 + rel_tol) );
    //         }

    //         if( (y_hot_L > 1 + rel_tol)  || (y_hot_L < 0))
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("y_hot_L not within (0,1 + rel_tol) \n");
    //             printf("y_hot_L = %e \n", y_hot_L);
    //             printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
    //             printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
    //             printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

    //             assert( (y_hot_L > 0) && (y_hot_L < 1 + rel_tol) );
    //         }

    //         if( (x_cold_L > 1 + rel_tol)  || (x_cold_L < 0))
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("x_cold_L not within (0,1 + rel_tol) \n");
    //             printf("x_cold_L = %e \n", x_cold_L);
    //             printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
    //             printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
    //             printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

    //             assert( (x_cold_L > 0) && (x_cold_L < 1 + rel_tol) );
    //         }

    //         if( (y_cold_L > 1 + rel_tol)  || (y_cold_L < 0))
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("y_cold_L not within (0,1 + rel_tol) \n");
    //             printf("y_cold_L = %e \n", y_cold_L);
    //             printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
    //             printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
    //             printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

    //             assert( (y_cold_L > 0) && (y_cold_L < 1 + rel_tol) );
    //         }


    //         if( std::abs(x_cold_L + x_hot_L - 1) >= rel_tol)
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("x_cold_L + x_hot_L =/= 1 \n");
    //             printf("x_cold_L           = %e \n", x_cold_L);
    //             printf("x_hot_L            = %e \n", x_hot_L);
    //             printf("x_cold_L + x_hot_L = %e \n", x_cold_L + x_hot_L);

    //             assert( std::abs(x_cold_L + x_hot_L - 1) < rel_tol);
    //         }

    //         if( std::abs(y_cold_L + y_hot_L - 1) >= rel_tol)
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("y_cold_L + y_hot_L =/= 1 \n");
    //             printf("y_cold_L           = %e \n", y_cold_L);
    //             printf("y_hot_L            = %e \n", y_hot_L);
    //             printf("y_cold_L + y_hot_L = %e \n", y_cold_L + y_hot_L);

    //             assert( std::abs(y_cold_L + y_hot_L - 1) < rel_tol);
    //         }

    // }
    // if (cR->multiphase)
    // {
    //     x_hot_R  = cR->cons_hot[DDD] / cR->cons[DDD];
    //     y_hot_R  = E_int_from_cons(cR->cons_hot) / E_int_from_cons(cR->cons);

    //     x_cold_R = cR->cons_cold[DDD] / cR->cons[DDD];
    //     y_cold_R = E_int_from_cons(cR->cons_cold) / E_int_from_cons(cR->cons);

    //         if( (x_hot_R > 1 + rel_tol)  || (x_hot_R < 0))
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("x_hot_R not within (0,1 + rel_tol) \n");
    //             printf("x_hot_R = %e \n", x_hot_R);
    //             printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
    //             printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
    //             printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

    //             assert( (x_hot_R > 0) && (x_hot_R < 1 + rel_tol) );
    //         }

    //         if( (y_hot_R > 1 + rel_tol)  || (y_hot_R < 0))
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("y_hot_R not within (0,1 + rel_tol) \n");
    //             printf("y_hot_R = %e \n", y_hot_R);
    //             printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
    //             printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
    //             printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

    //             assert( (y_hot_R > 0) && (y_hot_R < 1 + rel_tol) );
    //         }

    //         if( (x_cold_R > 1 + rel_tol)  || (x_cold_R < 0))
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("x_cold_R not within (0,1 + rel_tol) \n");
    //             printf("x_cold_R = %e \n", x_cold_R);
    //             printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
    //             printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
    //             printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

    //             assert( (x_cold_R > 0) && (x_cold_R < 1 + rel_tol) );
    //         }

    //         if( (y_cold_R > 1 + rel_tol)  || (y_cold_R < 0))
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("y_cold_R not within (0,1 + rel_tol) \n");
    //             printf("y_cold_R = %e \n", y_cold_R);
    //             printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
    //             printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
    //             printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

    //             assert( (y_cold_R > 0) && (y_cold_R < 1 + rel_tol) );
    //         }
            

    //         if( std::abs(x_cold_R + x_hot_R - 1) >= rel_tol)
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("x_cold_R + x_hot_R =/= 1 \n");
    //             printf("x_cold_R           = %e \n", x_cold_R);
    //             printf("x_hot_R            = %e \n", x_hot_R);
    //             printf("x_cold_R + x_hot_R = %e \n", x_cold_R + x_hot_R);

    //             assert( std::abs(x_cold_R + x_hot_R - 1) < rel_tol);
    //         }

    //         if( std::abs(y_cold_R + y_hot_R - 1) >= rel_tol)
    //         {
    //             printf("------ERROR in conduction()------- \n");
    //             printf("y_cold_R + y_hot_R =/= 1 \n");
    //             printf("y_cold_R           = %e \n", y_cold_R);
    //             printf("y_hot_R            = %e \n", y_hot_R);
    //             printf("y_cold_R + y_hot_R = %e \n", y_cold_R + y_hot_R);

    //             assert( std::abs(y_cold_R + y_hot_R - 1) < rel_tol);
    //         }
    // }

    // const double E_int_L_before = E_int_from_cons(cL->cons);
    // const double E_int_R_before = E_int_from_cons(cR->cons);

    // Apply conductive fluxes

    cL->cons[DDD] -= sgn * dM;
    cR->cons[DDD] += sgn * dM;

    cL->cons[SRR] -= sgn * dP;
    cR->cons[SRR] += sgn * dP;

    cL->cons[ZZZ] -= sgn * dM_Z;
    cR->cons[ZZZ] += sgn * dM_Z;

    cL->cons[TAU] -= sgn * dE; // dE already contains both kinetic and thermal energy
    cR->cons[TAU] += sgn * dE; // dE already contains both kinetic and thermal energy


    // const double E_int_L_after = E_int_from_cons(cL->cons);
    // const double E_int_R_after = E_int_from_cons(cR->cons);


    // question: should I only be applying fluxes from 
    // cold subgrid to hot subgrid, if possible?
    // That seems like it could lead to errors since the fluxes are 
    // computed using the total cell values
    // 
    // For now, transfer from/to both subgrid components
    if (cL->multiphase)
    {
        // const double dM_L = - sgn * dM;
        // const double dP_L = - sgn * dP;
        // const double dE_L = - sgn * dE;

        // const double dE_int_L = E_int_L_after - E_int_L_before;
        // const double dE_kin_L = dE_L - dE_int_L;

        // const double z_hot_L  = cL->cons_hot[ZZZ]  / cL->cons_hot[DDD];
        // const double z_cold_L = cL->cons_cold[ZZZ] / cL->cons_cold[DDD];

        cL->cons_hot[DDD]  -= cL->x_hot  * sgn * dM;      
        cL->cons_cold[DDD] -= cL->x_cold * sgn * dM;     

        cL->cons_hot[SRR]  -= cL->x_hot  * sgn * dP;      
        cL->cons_cold[SRR] -= cL->x_cold * sgn * dP;    

        cL->cons_hot[ZZZ]  -= cL->x_hot  * cL->z_hot  * sgn * dM;      
        cL->cons_cold[ZZZ] -= cL->x_cold * cL->z_cold * sgn * dM; 

        // cL->cons_hot[TAU]  += (y_hot_L  * dE_int_L) + (x_hot_L  * dE_kin_L);      
        // cL->cons_cold[TAU] += (y_cold_L * dE_int_L) + (x_cold_L * dE_kin_L); 

        // cL->cons_hot[TAU]  -= y_hot_L  * sgn * dE;      
        // cL->cons_cold[TAU] -= y_cold_L * sgn * dE; 
    }


    if (cR->multiphase)
    {

        // const double dM_R = + sgn * dM;
        // const double dP_R = + sgn * dP;
        // const double dE_R = + sgn * dE;

        // const double dE_int_R = E_int_R_after - E_int_R_before;
        // const double dE_kin_R = dE_R - dE_int_R;

        // const double z_hot_R  = cR->cons_hot[ZZZ]  / cR->cons_hot[DDD];
        // const double z_cold_R = cR->cons_cold[ZZZ] / cR->cons_cold[DDD];

        cR->cons_hot[DDD]  += cR->x_hot  * sgn * dM;      
        cR->cons_cold[DDD] += cR->x_cold * sgn * dM;      

        cR->cons_hot[SRR]  += cR->x_hot  * sgn * dP;      
        cR->cons_cold[SRR] += cR->x_cold * sgn * dP; 

        cR->cons_hot[ZZZ]  -= cR->x_hot  * cR->z_hot  * sgn * dM;      
        cR->cons_cold[ZZZ] -= cR->x_cold * cR->z_cold * sgn * dM; 

        // cR->cons_hot[TAU]  += (y_hot_R  * dE_int_R) + (x_hot_R  * dE_kin_R);      
        // cR->cons_cold[TAU] += (y_cold_R * dE_int_R) + (x_cold_R * dE_kin_R);

        // cR->cons_hot[TAU]  += y_hot_R  * sgn * dE;      
        // cR->cons_cold[TAU] += y_cold_R * sgn * dE; 
    }



    // ======== Verify post-conditions ========= //
    // #ifndef NDEBUG

    if (dM_dt_sat < dM_dt_unsat)
    {
        // not a post-condition, but a useful diagnostic
        printf("using saturated form of conduction\n");
        printf("T_L = %e \n", T_L);
        printf("T_R = %e \n", T_R);

        printf("dT^5/2 = %e \n", std::abs(std::pow(T_R, 5./2.) - std::pow(T_L, 5./2.)));

        printf("dM_dt_unsat = %e \n", dM_dt_unsat);
        printf("dM_dt_sat   = %e \n", dM_dt_sat);
    }
    // printf("sgn * x = %e \n", sgn * x);

    if (dM < 0)
    {
        printf("------ ERROR in conduction() postconditions------- \n");
        printf("Negative mass transfered! \n");
        printf("x = %e \n", x);
        assert(dM > 0);
    }

    if( x < 0)
    {
        printf("------ ERROR in conduction() postconditions------- \n");
        printf("Negative mass fraction transfered! \n");
        printf("x = %e \n", x);
        printf("T_L = %e \n", T_L);
        printf("T_R = %e \n", T_R);

        printf("dM_dt_unsat = %e \n", dM_dt_unsat);
        printf("dM_dt_sat   = %e \n", dM_dt_sat);
        assert(x >= 0);
    } 

    if( x >= 1)
    {
        printf("------ ERROR in conduction() postconditions------- \n");
        printf("Mass fraction >= 1 transfered! \n");
        printf("x = %e \n", x);
        printf("T_L = %e \n", T_L);
        printf("T_R = %e \n", T_R);

        printf("dM_dt_unsat = %e \n", dM_dt_unsat);
        printf("dM_dt_sat   = %e \n", dM_dt_sat);
        assert(x < 1);
    } 

    if (cL->multiphase)
    {
        if (cL->cons_hot[DDD] < 0)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[left] hot gas mass less than 0\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( cL->cons_hot[DDD] > 0 );
        }

        if (cL->cons_cold[DDD] < 0)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[left] cold gas mass less than 0\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( cL->cons_cold[DDD] > 0 );
        }

        if ( (cL->cons_hot[DDD] / cL->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[left] hot gas mass greater than total mass\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( (cL->cons_hot[DDD] / cL->cons[DDD]) < (1+rel_tol) );
        }

        if ( (cL->cons_cold[DDD] / cL->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[left] cold gas mass greater than total mass\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( (cL->cons_cold[DDD] / cL->cons[DDD]) < (1+rel_tol) );
        }

        if ( std::abs(1-( (cL->cons_cold[DDD] + cL->cons_hot[DDD])/cL->cons[DDD])) > rel_tol)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[left] cold mass + hot mass =/= total mass\n");
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert(  std::abs(1-( (cL->cons_cold[DDD] + cL->cons_hot[DDD])/cL->cons[DDD])) <= rel_tol);
        }




        if (cL->cons_hot[TAU] < 0)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[left] hot gas energy less than 0\n");
            printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

            assert( cL->cons_hot[TAU] > 0 );
        }

        if (cL->cons_cold[TAU] < 0)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[left] cold gas energy less than 0\n");
            printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

            assert( cL->cons_cold[TAU] > 0 );
        }

        // if (E_int_from_cons(cL->cons_hot) < 0)
        // {
        //     printf("------ERROR in conduction() postconditions------- \n");
        //     printf("[left] hot gas *internal* energy less than 0\n");
        //     printf("[left] E_int (hot)  = %e \n", E_int_from_cons(cL->cons_hot));
        //     printf("cL->cons_hot[TAU]   = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons[TAU]       = %e \n", cL->cons[TAU]);

        //     assert( E_int_from_cons(cL->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(cL->cons_cold) < 0)
        // {
        //     printf("------ERROR in conduction() postconditions------- \n");
        //     printf("[left] cold gas *internal* energy less than 0\n");
        //     printf("[left] E_int (cold)  = %e \n", E_int_from_cons(cL->cons_cold));
        //     printf("cL->cons_cold[TAU]   = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons[TAU]        = %e \n", cL->cons[TAU]);

        //     assert( E_int_from_cons(cL->cons_cold) > 0 );
        // }

        // if ( (cL->cons_hot[TAU] / cL->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in conduction() postconditions------- \n");
        //     printf("[left] hot gas energy greater than total energy\n");
        //     printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

        //     assert( (cL->cons_hot[TAU] / cL->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( (cL->cons_cold[TAU] / cL->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in conduction() postconditions------- \n");
        //     printf("[left] cold gas energy greater than total energy\n");
        //     printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

        //     assert( (cL->cons_cold[TAU] / cL->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( std::abs(1-( (cL->cons_cold[TAU] + cL->cons_hot[TAU])/cL->cons[TAU])) > rel_tol)
        // {
        //     printf("------ERROR in conduction() postconditions------- \n");
        //     printf("[left] cold energy + hot energy =/= total energy\n");
        //     printf("cL->cons_cold[TAU]                      = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons_hot[TAU]                       = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons_cold[TAU] + cL->cons_hot[TAU]  = %e \n", cL->cons_cold[TAU] + cL->cons_hot[TAU]);
        //     printf("cL->cons[TAU]                           = %e \n", cL->cons[TAU]);
        //     printf("relative error  = %e \n", 1 - ( (cL->cons_cold[TAU] + cL->cons_hot[TAU]) / cL->cons[TAU]));

        //     assert(  std::abs(1-( (cL->cons_cold[TAU] + cL->cons_hot[TAU])/cL->cons[TAU])) <= rel_tol);
        // }
    }

    if (cR->multiphase)
    {
        if (cR->cons_hot[DDD] < 0)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[right] hot gas mass less than 0\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( cR->cons_hot[DDD] > 0 );
        }

        if (cR->cons_cold[DDD] < 0)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[right] cold gas mass less than 0\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( cR->cons_cold[DDD] > 0 );
        }

        if ( (cR->cons_hot[DDD] / cR->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[right] hot gas mass greater than total mass\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( (cR->cons_hot[DDD] / cR->cons[DDD]) < (1+rel_tol) );
        }

        if ( (cR->cons_cold[DDD] / cR->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[right] cold gas mass greater than total mass\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( (cR->cons_cold[DDD] / cR->cons[DDD]) < (1+rel_tol) );
        }

        if ( std::abs(1-( (cR->cons_cold[DDD] + cR->cons_hot[DDD])/cR->cons[DDD])) > rel_tol)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[right] cold mass + hot mass =/= total mass\n");
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert(  std::abs(1-( (cR->cons_cold[DDD] + cR->cons_hot[DDD])/cR->cons[DDD])) <= rel_tol);
        }




        if (cR->cons_hot[TAU] < 0)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[right] hot gas energy less than 0\n");
            printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
            printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
            printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

            assert( cR->cons_hot[TAU] > 0 );
        }

        if (cR->cons_cold[TAU] < 0)
        {
            printf("------ERROR in conduction() postconditions------- \n");
            printf("[right] cold gas energy less than 0\n");
            printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
            printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
            printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

            assert( cR->cons_cold[TAU] > 0 );
        }

        // if (E_int_from_cons(cR->cons_hot) < 0)
        // {
        //     printf("------ERROR in riemann() postconditions------- \n");
        //     printf("[right] hot gas *internal* energy less than 0\n");
        //     printf("[right] E_int (hot) = %e \n", E_int_from_cons(cR->cons_hot));
        //     printf("cR->cons_hot[TAU]   = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons[TAU]       = %e \n", cR->cons[TAU]);

        //     assert( E_int_from_cons(cR->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(cR->cons_cold) < 0)
        // {
        //     printf("------ERROR in riemann() postconditions------- \n");
        //     printf("[right] cold gas *internal* energy less than 0\n");
        //     printf("[right] E_int (cold) = %e \n", E_int_from_cons(cR->cons_cold));
        //     printf("cR->cons_cold[TAU]   = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons[TAU]        = %e \n", cR->cons[TAU]);

        //     assert( E_int_from_cons(cR->cons_cold) > 0 );
        // }

        // if ( (cR->cons_hot[TAU] / cR->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in conduction() postconditions------- \n");
        //     printf("[right] hot gas energy greater than total energy\n");
        //     printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

        //     assert( (cR->cons_hot[TAU] / cR->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( (cR->cons_cold[TAU] / cR->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in conduction() postconditions------- \n");
        //     printf("[right] cold gas energy greater than total energy\n");
        //     printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

        //     assert( (cR->cons_cold[TAU] / cR->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( std::abs(1-( (cR->cons_cold[TAU] + cR->cons_hot[TAU])/cR->cons[TAU])) > rel_tol)
        // {
        //     printf("------ERROR in conduction() postconditions------- \n");
        //     printf("[right] cold energy + hot energy =/= total energy\n");
        //     printf("cR->cons_cold[TAU]                      = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons_hot[TAU]                       = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons_cold[TAU] + cR->cons_hot[TAU]  = %e \n", cR->cons_cold[TAU] + cR->cons_hot[TAU]);
        //     printf("cR->cons[TAU]                           = %e \n", cR->cons[TAU]);
        //     printf("relative error  = %e \n", 1 - ( (cR->cons_cold[TAU] + cR->cons_hot[TAU]) / cR->cons[TAU]));

        //     assert(  std::abs(1-( (cR->cons_cold[TAU] + cR->cons_hot[TAU])/cR->cons[TAU])) <= rel_tol);
        // }
    }


}



void subgrid_conduction( struct cell * c , const double dt )
{

    // ============================================= //
    //
    //  Calculate and add heat fluxes from subgrid conduction,
    //  using the equations of Keller et al. (2014), section 2.3 (eq. 12)
    //  This routine also adds in a saturation limit, which is not included
    //  in the original paper.
    //
    //  Inputs:
    //    - c         - cell
    //    - dt        - timestep [s]
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - overwrites cons_hot[DDD], cons_cold[DDD] (mass) 
    //      and cons_hot[TAU], cons_cold[TAU] (energy) 
    //
    //  Notes:
    //    - if cell not flagged as multiphase, nothing is changed
    //    - 
    //
    // ============================================= // 

    if (c->multiphase==0)
    {
        return;
    }

    // ======== Verify pre-conditions ========= //
    const double rel_tol = 1e-5; // relative tolerance for float comparison
    if (c->multiphase)
    {
        if (c->cons_hot[DDD] < 0)
        {
            printf("------ERROR in subgrid_conduction() preconditions------- \n");
            printf("hot gas mass less than 0\n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert( c->cons_hot[DDD] > 0 );
        }

        if (c->cons_cold[DDD] < 0)
        {
            printf("------ERROR in subgrid_conduction() preconditions------- \n");
            printf("cold gas mass less than 0\n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert( c->cons_cold[DDD] > 0 );
        }

        if ( (c->cons_hot[DDD] / c->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in subgrid_conduction() preconditions------- \n");
            printf("hot gas mass greater than total mass\n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert( (c->cons_hot[DDD] / c->cons[DDD]) < (1+rel_tol) );
        }

        if ( (c->cons_cold[DDD] / c->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in subgrid_conduction() preconditions------- \n");
            printf("cold gas mass greater than total mass\n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert( (c->cons_cold[DDD] / c->cons[DDD]) < (1+rel_tol) );
        }

        if ( std::abs(1-( (c->cons_cold[DDD] + c->cons_hot[DDD])/c->cons[DDD])) > rel_tol)
        {
            printf("------ERROR in subgrid_conduction() preconditions------- \n");
            printf("cold mass + hot mass =/= total mass\n");
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert(  std::abs(1-( (c->cons_cold[DDD] + c->cons_hot[DDD])/c->cons[DDD])) <= rel_tol);
        }




        if (c->cons_hot[TAU] < 0)
        {
            printf("------ERROR in subgrid_conduction() preconditions------- \n");
            printf("hot gas energy less than 0\n");
            printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
            printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
            printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

            assert( c->cons_hot[TAU] > 0 );
        }

        if (c->cons_cold[TAU] < 0)
        {
            printf("------ERROR in subgrid_conduction() preconditions------- \n");
            printf("cold gas energy less than 0\n");
            printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
            printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
            printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

            assert( c->cons_cold[TAU] > 0 );
        }

        // if (E_int_from_cons(c->cons_hot) < 0)
        // {
        //     printf("------ERROR in subgrid_conduction() preconditions------- \n");
        //     printf("hot gas *internal* energy less than 0\n");
        //     printf("E_int (hot)  = %e \n", E_int_from_cons(c->cons_hot));
        //     printf("c->cons_hot[TAU]   = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons[TAU]       = %e \n", c->cons[TAU]);

        //     assert( E_int_from_cons(c->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(c->cons_cold) < 0)
        // {
        //     printf("------ERROR in subgrid_conduction() preconditions------- \n");
        //     printf("cold gas *internal* energy less than 0\n");
        //     printf("E_int (cold)  = %e \n", E_int_from_cons(c->cons_cold));
        //     printf("c->cons_cold[TAU]   = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]        = %e \n", c->cons[TAU]);

        //     assert( E_int_from_cons(c->cons_cold) > 0 );
        // }

        // if ( (c->cons_hot[TAU] / c->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in subgrid_conduction() preconditions------- \n");
        //     printf("hot gas energy greater than total energy\n");
        //     printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

        //     assert( (c->cons_hot[TAU] / c->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( (c->cons_cold[TAU] / c->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in subgrid_conduction() preconditions------- \n");
        //     printf("cold gas energy greater than total energy\n");
        //     printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

        //     assert( (c->cons_cold[TAU] / c->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( std::abs(1-( (c->cons_cold[TAU] + c->cons_hot[TAU])/c->cons[TAU])) > rel_tol)
        // {
        //     printf("------ERROR in subgrid_conduction() preconditions------- \n");
        //     printf("cold energy + hot energy =/= total energy\n");
        //     printf("c->cons_cold[TAU]                      = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons_hot[TAU]                       = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons_cold[TAU] + c->cons_hot[TAU]  = %e \n", c->cons_cold[TAU] + c->cons_hot[TAU]);
        //     printf("c->cons[TAU]                           = %e \n", c->cons[TAU]);
        //     printf("relative error  = %e \n", 1 - ( (c->cons_cold[TAU] + c->cons_hot[TAU]) / c->cons[TAU]));

        //     assert(  std::abs(1-( (c->cons_cold[TAU] + c->cons_hot[TAU])/c->cons[TAU])) <= rel_tol);
        // }
    }


    // ======== Primary Code ========= //


    const double T_hot = calc_T(c->prim_hot);

    const double dr = c->dr;
    const double mu = get_mean_molecular_weight(c->prim[ZZZ]);


    // calculate *UNsaturated* flux (the original eq 12)
    const double dM_b_dt_unsat = 16.*M_PI*mu*m_proton / (25.*k_boltzmann)
        * kappa_0 * std::pow(T_hot, 5./2.) * dr;




    // calculate *saturated* flux (my addition to eq 12)
    // const double V_cold = c->cons_cold[DDD] / c->cons_cold[RHO];
    const double A_cold = 4 * std::pow((3. * c->V_cold / (4. * M_PI)), 2./3.);

    const double phi_s = 1.1; // see note near Eq. 8 of Cowie & McKee 1977
    const double c_s_cold = std::sqrt(std::abs(GAMMA_LAW*c->prim_cold[PPP]/c->prim_cold[RHO]));

    const double q = 5 * phi_s * c_s_cold * c->prim_cold[PPP];
    const double dM_b_dt_sat = q * (A_cold) / (c_s_cold*c_s_cold); // it's the colder cell's energy which is transported


    // choose saturated or unsaturated
    const double dM_b_dt = std::min(dM_b_dt_unsat, dM_b_dt_sat);

    const double dM_b = dM_b_dt * dt;


    // to do: also transfer the correct amount of kinetic energy!!!

    const double x = dM_b / c->cons_cold[DDD];
    const double dE = dM_b * c_s_cold * c_s_cold;
    const double dM_Z = x * c->cons_cold[ZZZ];

    // ======== Verify before transferring mass / energy ========= //
    // #ifndef NDEBUG

    if (!std::isfinite(dM_b_dt_unsat))
    {
        printf("------ ERROR in subgrid_conduction()------- \n");
        printf("Non-finite unsaturated conduction rate! \n");
        printf("dM_b_dt_unsat = %e \n", dM_b_dt_unsat);
        printf("T_hot         = %e \n", T_hot);
        printf("T_hot**(5/2)  = %e \n", std::pow(T_hot, 5./2.));
        printf("dr            = %e \n", dr);
        assert( std::isfinite(dM_b_dt_unsat) );

    }

    if (!std::isfinite(dM_b_dt_sat))
    {
        printf("------ ERROR in subgrid_conduction()------- \n");
        printf("Non-finite saturated conduction rate! \n");
        printf("dM_b_dt_sat       = %e \n", dM_b_dt_sat);
        printf("c->V_cold         = %e \n", c->V_cold);
        printf("c->V_hot          = %e \n", c->V_hot);
        printf("A_cold            = %e \n", A_cold);
        printf("c->prim_cold[PPP] = %e \n", c->prim_cold[PPP]);
        printf("c->prim_cold[RHO] = %e \n", c->prim_cold[RHO]);
        printf("c_s_cold          = %e \n", c_s_cold);
        printf("q                 = %e \n", q);
        assert( std::isfinite(dM_b_dt_sat) );

    }

    if (dM_b < 0)
    {
        printf("------ ERROR in subgrid_conduction()------- \n");
        printf("Negative mass transfered! \n");
        printf("dM_b = %e \n", dM_b);
        assert(dM_b > 0);
    }

    if( x < 0)
    {
        printf("------ ERROR in subgrid_conduction()------- \n");
        printf("Negative mass fraction transfered! \n");
        printf("x = %e \n", x);
        printf("T_hot = %e \n", T_hot);

        printf("dM_b_dt = %e \n", dM_b_dt);
        printf("dt      = %e \n", dt);
        printf("dM_b    = %e \n", dM_b);

        assert(x >= 0);
    } 

    if( x >= 1)
    {
        const double T_cold = calc_T(c->prim_cold);
        printf("------ ERROR in subgrid_conduction()------- \n");
        printf("Mass fraction >= 1 transfered! \n");
        printf("x       = %e \n", x);
        printf("T_hot   = %e \n", T_hot);
        printf("T_cold  = %e \n", T_cold);
        printf("dr      = %e \n", dr);
        printf("mu      = %e \n", mu);


        printf("dM_b_dt_sat   = %e \n", dM_b_dt_sat);
        printf("dM_b_dt_unsat = %e \n", dM_b_dt_unsat);
        printf("dM_b_dt = %e \n", dM_b_dt);
        printf("dt      = %e \n", dt);
        printf("dM_b    = %e \n", dM_b);
        printf("M_hot   = %e \n", c->cons_hot[DDD]);
        printf("M_cold  = %e \n", c->cons_cold[DDD]);


        assert(x < 1);
    } 

    if (c->cons_cold[TAU] <= 0)
    {

        printf("------ ERROR in subgrid_conduction()------- \n");
        printf("No cold energy to begin with! \n");
        printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        assert(0 < c->cons_cold[TAU]);
    }   


    if (dE >= c->cons_cold[TAU])
    {
        const double T_cold = calc_T(c->prim_cold);

        printf("------ ERROR in subgrid_conduction()------- \n");
        printf("All energy transfered! \n");
        printf("dE                = %e \n", dE);
        printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        printf("T_hot   = %e \n", T_hot);
        printf("T_cold  = %e \n", T_cold);

        assert(dE < c->cons_cold[TAU]);
    }   




    // Apply conductive fluxes

    c->cons_cold[DDD] -= dM_b;
    c->cons_hot[DDD]  += dM_b;

    c->cons_cold[ZZZ] -= dM_Z;
    c->cons_hot[ZZZ]  += dM_Z;

    c->cons_cold[TAU] -= dE;
    c->cons_hot[TAU]  += dE;


    // ======== Verify post-conditions ========= //
    if (c->multiphase)
    {
        if (c->cons_hot[DDD] < 0)
        {
            printf("------ERROR in subgrid_conduction() postconditions------- \n");
            printf("hot gas mass less than 0\n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert( c->cons_hot[DDD] > 0 );
        }

        if (c->cons_cold[DDD] < 0)
        {
            printf("------ERROR in subgrid_conduction() postconditions------- \n");
            printf("cold gas mass less than 0\n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert( c->cons_cold[DDD] > 0 );
        }

        if ( (c->cons_hot[DDD] / c->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in subgrid_conduction() postconditions------- \n");
            printf("hot gas mass greater than total mass\n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert( (c->cons_hot[DDD] / c->cons[DDD]) < (1+rel_tol) );
        }

        if ( (c->cons_cold[DDD] / c->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in subgrid_conduction() postconditions------- \n");
            printf("cold gas mass greater than total mass\n");
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert( (c->cons_cold[DDD] / c->cons[DDD]) < (1+rel_tol) );
        }

        if ( std::abs(1-( (c->cons_cold[DDD] + c->cons_hot[DDD])/c->cons[DDD])) > rel_tol)
        {
            printf("------ERROR in subgrid_conduction() postconditions------- \n");
            printf("cold mass + hot mass =/= total mass\n");
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);

            assert(  std::abs(1-( (c->cons_cold[DDD] + c->cons_hot[DDD])/c->cons[DDD])) <= rel_tol);
        }




        if (c->cons_hot[TAU] < 0)
        {
            printf("------ERROR in subgrid_conduction() postconditions------- \n");
            printf("hot gas energy less than 0\n");
            printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
            printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
            printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

            assert( c->cons_hot[TAU] > 0 );
        }

        if (c->cons_cold[TAU] < 0)
        {
            printf("------ERROR in subgrid_conduction() postconditions------- \n");
            printf("cold gas energy less than 0\n");
            printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
            printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
            printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

            assert( c->cons_cold[TAU] > 0 );
        }

        // if (E_int_from_cons(c->cons_hot) < 0)
        // {
        //     printf("------ERROR in subgrid_conduction() postconditions------- \n");
        //     printf("hot gas *internal* energy less than 0\n");
        //     printf("E_int (hot)  = %e \n", E_int_from_cons(c->cons_hot));
        //     printf("c->cons_hot[TAU]   = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons[TAU]       = %e \n", c->cons[TAU]);

        //     assert( E_int_from_cons(c->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(c->cons_cold) < 0)
        // {
        //     printf("------ERROR in subgrid_conduction() postconditions------- \n");
        //     printf("cold gas *internal* energy less than 0\n");
        //     printf("E_int (cold)  = %e \n", E_int_from_cons(c->cons_cold));
        //     printf("c->cons_cold[TAU]   = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]        = %e \n", c->cons[TAU]);

        //     assert( E_int_from_cons(c->cons_cold) > 0 );
        // }

        // if ( (c->cons_hot[TAU] / c->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in subgrid_conduction() postconditions------- \n");
        //     printf("hot gas energy greater than total energy\n");
        //     printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

        //     assert( (c->cons_hot[TAU] / c->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( (c->cons_cold[TAU] / c->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in subgrid_conduction() postconditions------- \n");
        //     printf("cold gas energy greater than total energy\n");
        //     printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

        //     assert( (c->cons_cold[TAU] / c->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( std::abs(1-( (c->cons_cold[TAU] + c->cons_hot[TAU])/c->cons[TAU])) > rel_tol)
        // {
        //     printf("------ERROR in subgrid_conduction() postconditions------- \n");
        //     printf("cold energy + hot energy =/= total energy\n");
        //     printf("c->cons_cold[TAU]                      = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons_hot[TAU]                       = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons_cold[TAU] + c->cons_hot[TAU]  = %e \n", c->cons_cold[TAU] + c->cons_hot[TAU]);
        //     printf("c->cons[TAU]                           = %e \n", c->cons[TAU]);
        //     printf("relative error  = %e \n", 1 - ( (c->cons_cold[TAU] + c->cons_hot[TAU]) / c->cons[TAU]));

        //     assert(  std::abs(1-( (c->cons_cold[TAU] + c->cons_hot[TAU])/c->cons[TAU])) <= rel_tol);
        // }
    }



}

void vel( const double * prim1 , const double * prim2 , 
          double * Sl , double * Sr , double * Ss )
{

    // ============================================= //
    //
    //  Calculate linearized wave family velocities
    //
    //  Inputs:
    //    - prim1     - primitive variable array for INNER cell
    //    - prim2     - primitive variable array for OUTER cell
    //    - Sl        -  left-most moving sound wave
    //    - Sr        - right-most moving sound wave
    //    - Ss        - entropy wave
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - overwrites Sl, Sr, Ss as output variables
    //
    //  Notes:
    //    - Nomenclature assumes "left" is "inner" (smaller radius)
    //    - The physics doesn't care which is prim1 or prim2,
    //      but it does matter for the sign convention on velocities
    //
    // ============================================= // 

    const double gam  = GAMMA_LAW;

    const double P1   = prim1[PPP];
    const double rho1 = prim1[RHO];
    const double vn1  = prim1[VRR];

    const double cs1  = std::sqrt(std::abs(gam*P1/rho1));

    const double P2   = prim2[PPP];
    const double rho2 = prim2[RHO];
    const double vn2  = prim2[VRR];

    const double cs2  = std::sqrt(std::abs(gam*P2/rho2));

    *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

    *Sr =  cs1 + vn1;
    *Sl = -cs1 + vn1;

    if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
    if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;

    // ======== Verify post-conditions ========= //
    // wave family characteristics should have different speeds
    #ifndef NDEBUG
    if( *Sr == *Sl)
    {
        printf("------- ERROR in vel() ----------- \n");
        printf("Sr = Sl = %e \n", *Sr);
        printf("cs1  = %e \n", cs1);
        printf("cs2  = %e \n", cs2);
        printf("vn1  = %e \n", vn1);
        printf("vn1  = %e \n", vn2);
        printf("P1   = %e \n", P1);
        printf("rho1 = %e \n", rho1);
        printf("P2   = %e \n", P2);
        printf("rho2 = %e \n", rho2);
        assert( *Sr != *Sl );
    }
    #endif
  
}

double mindt( const double * prim , const double w , 
              const double r , const double dr )
{

    // ============================================= //
    //
    //  For a single cell, find the minimum wave crossing time,
    //  in the cell's rest-frame
    //
    //  Inputs:
    //    - prim      - primitive variable array for INNER cell
    //    - w         - cell boundary velocity, to subtract from wave families
    //    - r         - radius of cell center
    //                - only used if Rayleigh-Taylor code is active
    //    - dr        - width of cell
    //
    //  Returns:
    //       dt       - minimum wave-crossing time of the 3 waves
    //
    //  Side effects:
    //    None
    //
    //  Notes:
    //    - w is probably the AVERAGE of the velocity of both boundaries
    //      (see getmindt())
    //
    // ============================================= // 

    const double rho = prim[RHO];
    const double P   = prim[PPP];
    const double vr  = prim[VRR];
    const double gam = GAMMA_LAW;

    const double cs  = std::sqrt(std::abs(gam*P/rho));

    const double maxvr  = cs + std::abs( vr - w );
    const double dt     = dr/maxvr;

    return dt;
}


