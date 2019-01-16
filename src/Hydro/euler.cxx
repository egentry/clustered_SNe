
#include <cmath>
#include <assert.h>
#include <stdlib.h>
#include <functional>
extern "C" {
#include <grackle.h>
}
#include <lapacke.h>

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
                const double dV , const bool verbose
                , const std::string called_from
                )
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
            printf("called from: %s\n", called_from.c_str());
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
    
    verify_multiphase_conditions(c, "calc_multiphase_prim()", "pre");

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

double E_kin_from_cons( const double * cons )
{

    const double E_kin = .5*cons[SRR]*cons[SRR] / cons[DDD];

    return E_kin;
}

double E_int_from_cons( const double * cons )
{

    // const double mass  = cons[DDD];
    // const double vr    = cons[SRR] / mass;
    const double E_int = cons[TAU] - E_kin_from_cons(cons);

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

        c->cons_hot[SRR]  = c->cons[SRR] * (c->cons_hot[DDD]  / c->cons[DDD]);
        c->cons_cold[SRR] = c->cons[SRR] * (c->cons_cold[DDD] / c->cons[DDD]);

        // assert(E_int_from_cons(c->cons_cold) > 0);
        // assert(E_int_from_cons(c->cons_hot) > 0);
        // assert(E_int_from_cons(c->cons) > 0);
    }       

    // now decompose the change in energy, and apply it to subgrid components
    if (c->multiphase)
    {
        const double E_int = E_int_from_cons(c->cons);
        const double E_kin = E_kin_from_cons(c->cons);

        const double dE_int = E_int - c->E_int_initial;
        const double dE_kin = E_kin - c->E_kin_initial;

        printf("------- In source (*before* apply decomposed energy) ---------- \n");
        printf("E_int (hot)       = %e \n", E_int_from_cons(c->cons_hot));
        printf("E_int (cold)      = %e \n", E_int_from_cons(c->cons_cold));
        printf("E_int (total)     = %e \n", E_int_from_cons(c->cons));
        printf("c->cons_hot[SRR]  = %e \n", c->cons_hot[SRR]);
        printf("c->cons_cold[SRR] = %e \n", c->cons_cold[SRR]);

        printf("\n");
        printf("dE_int = %e\n", dE_int);
        printf("dE_kin = %e\n", dE_kin);
        printf("\n");

        // const double x_hot  = c->cons_hot[DDD] / c->cons[DDD];
        // const double x_cold = c->cons_cold[DDD] / c->cons[DDD];

        c->cons_hot[TAU]  += (c->y_hot  * dE_int) + (c->x_hot  * dE_kin);      
        c->cons_cold[TAU] += (c->y_cold * dE_int) + (c->x_cold * dE_kin);

        printf("------- In source (*after* apply decomposed energy) ---------- \n");
        printf("E_int (hot)       = %e \n", E_int_from_cons(c->cons_hot));
        printf("E_int (cold)      = %e \n", E_int_from_cons(c->cons_cold));
        printf("E_int (total)     = %e \n", E_int_from_cons(c->cons));
        printf("c->cons_hot[SRR]  = %e \n", c->cons_hot[SRR]);
        printf("c->cons_cold[SRR] = %e \n", c->cons_cold[SRR]);

        if(1)
        {

            printf("\n");


            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_hot[SRR]  = %e \n", c->cons_hot[SRR]);
            printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
            printf("c->cons_hot[ZZZ]  = %e \n", c->cons_hot[ZZZ]);

            printf("\n");

            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
            printf("c->cons_cold[SRR] = %e \n", c->cons_cold[SRR]);
            printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
            printf("c->cons_cold[ZZZ] = %e \n", c->cons_cold[ZZZ]);

            printf("\n");

            printf("c->cons[DDD]      = %e \n", c->cons[DDD]);
            printf("c->cons[SRR]      = %e \n", c->cons[SRR]);
            printf("c->cons[TAU]      = %e \n", c->cons[TAU]);
            printf("c->cons[ZZZ]      = %e \n", c->cons[ZZZ]);

            printf("\n");

            printf("E_kin (hot)       = %e \n", c->cons_hot[TAU]  - E_int_from_cons(c->cons_hot));
            printf("E_kin (cold)      = %e \n", c->cons_cold[TAU] - E_int_from_cons(c->cons_cold));
            printf("E_kin (total)     = %e \n", c->cons[TAU]      - E_int_from_cons(c->cons));
            printf("\n");

            printf("VRR (hot)         = %e \n", c->cons_hot[SRR]  / c->cons_hot[DDD]);
            printf("VRR (cold)        = %e \n", c->cons_cold[SRR] / c->cons_cold[DDD]);
            printf("VRR (total)       = %e \n", c->cons[SRR]      / c->cons[DDD]);
            printf("\n");
        }

        printf("before asserts\n");

        printf("E_int_from_cons(c->cons_hot) = %e \n", E_int_from_cons(c->cons_hot)); 
        printf("E_int_from_cons(c->cons_cold) = %e \n", E_int_from_cons(c->cons_cold));
        printf("E_int_from_cons(c->cons) = %e \n", E_int_from_cons(c->cons)); 

        assert(E_int_from_cons(c->cons_hot)  > 0); 
        assert(E_int_from_cons(c->cons_cold) > 0);
        assert(E_int_from_cons(c->cons)      > 0); 

        printf("past asserts\n");

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

            const bool allow_caching = false;
            const double dE_cool_hot  = cooling->calc_cooling(c->prim_hot,  c->cons_hot,  dt , cached_cooling , allow_caching ) * c->V_hot;
            const double dE_cool_cold = cooling->calc_cooling(c->prim_cold, c->cons_cold, dt , cached_cooling , allow_caching ) * c->V_cold;

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

            printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
            printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]); 
            printf("c->cons[DDD]      = %e \n", c->cons[DDD]); 

            printf("c->cons_hot[SRR]  = %e \n", c->cons_hot[SRR]);
            printf("c->cons_cold[SRR] = %e \n", c->cons_cold[SRR]); 
            printf("c->cons[SRR]      = %e \n", c->cons[SRR]); 
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
            bool allow_caching = true;
            double dr = rp - rm;
            if( r > (R_shock+(10*dr)) ) cached_cooling = true;

            c->dE_cool   = cooling->calc_cooling(c->prim, c->cons, dt , 
                                                 cached_cooling , 
                                                 allow_caching ) * dV;
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

void artificial_conduction( struct cell * cL , struct cell * cR, 
                            const double dA , const double dt )
{

    // ============================================= //
    //
    //  Calculate and add heat fluxes from artifical conduction,
    //  using the formalism of Noh (1987), esp. Eq 2.3
    //
    //  Inputs:
    //    - cL        - cell to the Left
    //    - prim2     - cell to the Right
    //    - dA        - interface area between cells [cm^2]
    //    - dt        - timestep [s]
    //
    //  Static Variables:
    //    - H_0 - shock conduction strength (scales like delta_u)
    //    - H_1 - thermal conduction strength (scales like sound speed) 
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - overwrites cons[TAU] (energy) of both cells cL, cR
    //      (cL, cR on left and right side of interface, respectively)
    //
    //  Notes:
    //    - The physics doesn't care which is cL or cR,
    //      but it *does* matter for the sign convention on velocities
    //      - and velocity sign convention sets if it's compressing
    //        (if not compressing, returns without adding flux)
    //
    // ============================================= // 

    double prim_L[NUM_Q];
    double prim_R[NUM_Q];

    const double dr_L = .5*cL->dr;
    const double dr_R = .5*cR->dr;

    for( int q=0 ; q<NUM_Q ; ++q )
    {
        prim_L[q] = cL->prim[q] + cL->grad[q]*dr_L;
        prim_R[q] = cR->prim[q] - cR->grad[q]*dr_R;
    }

    const double delta_u = prim_R[VRR] - prim_L[VRR];
    if (delta_u > 0) return;

    const double rho_average = 2 * (prim_L[RHO] * prim_R[RHO]) 
                                 / (prim_L[RHO] + prim_R[RHO]);

    // change in specific internal energy (E_int / mass)
    // it's missing a factor of (gamma-1), but that's just a constant term
    const double delta_E = ((prim_R[PPP]/prim_R[RHO]) - (prim_L[PPP]/prim_L[RHO])) / 2;

    const double cs_L = std::sqrt(std::abs(GAMMA_LAW*prim_L[PPP]/prim_L[RHO]));
    const double cs_R = std::sqrt(std::abs(GAMMA_LAW*prim_R[PPP]/prim_R[RHO]));
    const double cs_average = 2 * (cs_L * cs_R)
                                / (cs_L + cs_R);

    // minus sign between the two terms (which doesn't show up in Noh),
    // because of sign conventions + absolute values
    const double H_L = rho_average * ( (H_0 * delta_u) - (H_1 * cs_average) ) * delta_E;

    cL->cons[TAU] -= H_L * dA * dt;
    cR->cons[TAU] += H_L * dA * dt;

}

int index_2d_to_1d(const int i, const int j, const int I_max, const int J_max)
{
    const int index_1d = (i*J_max) + j;
    return index_1d;
}

void thermal_conduction_apply_linsolve( const double * e_int_old ,
                                const double * e_int_guess ,
                                double * e_int_new ,
                                const double * r ,
                                const double * V ,
                                const double * rho ,
                                const double * mu ,
                                const double dt ,
                                const int num_cells ,
                                const int verbose )
{

    // ============================================= //
    //
    //  Linearizes and solves the conduction equation for an implicit _guess_
    //  of e_int_guess. Returns a linearized solution e_int_new. 
    //  For Picard iteration, e_int_new becomes e_int_guess of the next iteration.
    //
    //  Inputs:
    //    - e_int_old, e_int_guess, e_int_new
    //        - array of internal energy per unit mass
    //        - has length num_cells
    //        - e_int_new is an output, and can be garbage on input
    //    - r  - array of cell-centered radii for each cell
    //    - V   - array of cell volumes
    //    - mu  - array of cell mean molecular weights
    //    - rho - array of cell densities
    //    - dt  - timestep
    //    - num_cells - the number of cells considered, not including guard cells
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - e_int_new completely overwritten with result
    //
    //  Notes:
    //    - The banded matrix solver we use is LAPACKE_dgbsv
    //      (This uses the c binding for the fortran LAPACK package)
    //    - For a good example of using LAPACK from c, and for the 
    //      indexing convention, check out:
    //      http://www.netlib.org/lapack/explore-html/d8/dd5/example___d_g_e_l_s__rowmajor_8c_source.html
    //    - I chose a banded diagonal solver when I was trying to solve a more
    //      complex problem. In principle, you could go back to using a simpler
    //      tri-diagonal solver, but I've got this debugged already
    //
    // ============================================= // 

    if(verbose)
    {
        printf("just entered thermal_conduction_apply_linsolve\n");
        fflush(stdout);
    }



    const int system_rank = num_cells;

    // Declare and initialize matrices so we can do:
    // A x = b
    // where b is `e_int_old`
    // and A is determined using `e_int_guess`

    double *A_matrix = new double[system_rank*system_rank];

    // initialize linear algebra containers
    for( int i=0 ; i<system_rank ; ++i )
    {
        for( int j=0 ; j<system_rank ; ++j )
        {
            A_matrix[index_2d_to_1d(i, j, system_rank, system_rank)] = 0.;
        }
    }

    if(verbose)
    {
        printf("adding identity matrix\n");
        fflush(stdout);
    }

    // add identity along diagonal
    for( int i=0 ; i<system_rank ; ++i )
    {
        A_matrix[index_2d_to_1d(i, i, system_rank, system_rank)] = 1.;
    }

    for( int i=0 ; i<(num_cells-1) ; ++i )
    {
        // this is interface i + 1/2
        const double r_interface = (r[i+1] + r[i]) / 2.;
        const double dA = get_dA(r_interface);
        const double dr = r[i+1] - r[i];

        const double mu_L = mu[i];
        const double mu_R = mu[i+1];

        const double e_int_L = e_int_guess[i];
        const double e_int_R = e_int_guess[i+1];

        const double T_L = (mu_L * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_L;
        const double T_R = (mu_R * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_R;

        // first, check if the temperatures are relatively close
        // then, check if there are effectively equal
        // if so, skip the conduction routine 
        //    (if you don't you might get silent divide-by-zero errors leading to NaN energy flux)
        if (std::abs(std::log10(T_L/T_R)) < 1)
        {
            if (std::abs(1 - (T_L/T_R)) < .01)
            {
                // effectively no temperature gradient; don't apply conduction
                // otherwise you get silent divide-by-zero errors, leading to 
                // fluxes that are NaN's
                continue;
            }
        }

        // okay, so now we need to know if it's going to be saturated or unsaturated.
        // The only real way to do this is actually calculate the energy flux.

        double rho_upwind;
        double e_int_upwind;

        int sgn; // sign of dT/dr deriv; (+1 if energy is transfered left, -1 if energy transfered right)
        int i_upwind;
        int i_downwind;
        if( e_int_L > e_int_R )
        {   
            i_upwind   = i;
            i_downwind = i+1;
            sgn        = -1;
            // left cell is hotter and therefore "upwind"
            rho_upwind = rho[i_upwind];
            e_int_upwind = e_int_guess[i_upwind];
        }
        else
        {
            i_upwind   = i + 1;
            i_downwind = i;
            sgn        = +1;
            // right cell is hotter and therefore "upwind"
            rho_upwind = rho[i_upwind];
            e_int_upwind = e_int_guess[i_upwind];
        }

        const double mu_average = (mu_L + mu_R) / 2.;
        const double e_int_average = (e_int_L + e_int_R) / 2.;

        const double Q_unsat = - kappa_0 * (2./7.) 
                                         * std::pow(mu_average * m_proton * (GAMMA_LAW-1) / k_boltzmann, 7./2.)
                                         * (std::pow(e_int_R, 7./2.) - std::pow(e_int_L, 7./2.) ) 
                                         / dr;




        // calculate *SATURATED* flux
        const double phi_s = 1.1; // see note near Eq. 8 of Cowie & McKee 1977
        const double Q_sat = - sgn * 5 
                                   * phi_s 
                                   * rho_upwind
                                   * std::pow( (GAMMA_LAW-1) * e_int_upwind, 3./2.) ;

        if( !std::isfinite(Q_sat) )
        {
            printf("Value Error: Q_sat =  %e for interface between cells %d and %d\n", 
                Q_sat, i, i+1);
            printf("i_upwind = %d\n", i_upwind);
            printf("sgn = %d\n", sgn);
            printf("phi_s = %e\n", phi_s);
            printf("rho_upwind= %e\n", rho_upwind);
            printf("e_int_upwind = %e\n", e_int_upwind);
            printf("T_L = %e\n", T_L);
            printf("T_R = %e\n", T_R);

            printf("\n");
            printf("e_int_guess[i] = %e\n", e_int_guess[i]);
            printf("e_int_guess[i+1] = %e\n", e_int_guess[i+1]);

            fflush(stdout);

            assert(std::isfinite(Q_sat));
        }

        if( !std::isfinite(Q_unsat) )
        {
            printf("Value Error: Q_unsat =  %e for interface between cells %d and %d\n", 
                Q_unsat, i, i+1);
            printf("i_upwind = %d\n", i_upwind);
            printf("sgn = %d\n", sgn);

            printf("T_L = %e\n", T_L);
            printf("T_R = %e\n", T_R);

            printf("std::pow(T_R, 7./2.) - std::pow(T_L, 7./2.) = %e\n", std::pow(T_R, 7./2.) - std::pow(T_L, 7./2.));            
            printf("dr = %e\n", dr);
            printf("e_int_upwind = %e\n", e_int_upwind);

            printf("\n");
            printf("rho[i]           = %e\n", rho[i]);
            printf("e_int_guess[i]   = %e\n", e_int_guess[i]);

            printf("rho[i+1]         = %e\n", rho[i+1]);
            printf("e_int_guess[i+1] = %e\n", e_int_guess[i+1]);

            fflush(stdout);

            assert(std::isfinite(Q_unsat));
        }


        if( std::abs(Q_sat) > std::abs(Q_unsat))
        {
            // UNSATURATED REGIME

            if( i == 0 )
            {
                if(verbose)
                {
                    printf("Interface between cells 0 + 1 unsaturated\n");
                    printf("Q_sat   = %e\n", Q_sat);
                    printf("Q_unsat = %e\n", Q_unsat);
                    printf("T_L = %e\n", T_L);
                    printf("T_R = %e\n", T_R);
                }
            }


            const double C_unsat = kappa_0 
                * std::pow( mu_average * m_proton * (GAMMA_LAW-1) / k_boltzmann , 7./2.)
                * std::pow(e_int_average, 5./2.) ;





            A_matrix[index_2d_to_1d(i, i, 
                                    system_rank, system_rank)] += C_unsat
                * dA
                * dt
                / (V[i] * rho[i])
                / dr;

            A_matrix[index_2d_to_1d(i, i+1, 
                                    system_rank, system_rank)] -= C_unsat
                * dA
                * dt
                / (V[i] * rho[i])
                / dr;


            A_matrix[index_2d_to_1d(i+1, i, 
                                    system_rank, system_rank)] -= C_unsat
                * dA
                * dt
                / (V[i+1] * rho[i+1])
                / dr;


            A_matrix[index_2d_to_1d(i+1, i+1, 
                                    system_rank, system_rank)] += C_unsat
                * dA
                * dt
                / (V[i+1] * rho[i+1])
                / dr;

        }
        else
        {
            // SATURATED REGIME
            // Okay, this one is easy. We don't have to recompute anything
            // because the matrix is just going to be F_mass * (U/rho)_upwind

            if( i == 0 )
            {
                if(verbose)
                {
                    printf("Interface between cells 0 + 1 saturated\n");
                    printf("Q_sat   = %e\n", Q_sat);
                    printf("Q_unsat = %e\n", Q_unsat);
                    printf("T_L = %e\n", T_L);
                    printf("T_R = %e\n", T_R);
                    printf("dE = %e\n", Q_sat * dA * dt );
                    printf("E_upwind   = %e\n", V[i_upwind]   * rho[i_upwind]   * e_int_old[i_upwind]);
                    printf("E_downwind = %e\n", V[i_downwind] * rho[i_downwind] * e_int_old[i_downwind]);

                    printf("\n");

                    printf("phi_s = %e\n", phi_s);
                    printf("rho[i_upwind]         = %e\n", rho[i_upwind]);
                    printf("rho[i_downwind]       = %e\n", rho[i_downwind]);
                    printf("e_int_old[i_upwind]   = %e\n", e_int_old[i_upwind]);
                    printf("e_int_old[i_downwind] = %e\n", e_int_old[i_downwind]);

                    printf("V_upwind = %e\n", V[i_upwind]);

                    printf("dt = %e\n", dt);
                    printf("dA = %e\n", dA);
                }
            }

            const double C_sat = 5 * phi_s
                                   * std::pow(GAMMA_LAW-1, 3./2.)
                                   * rho_upwind
                                   * std::pow(e_int_upwind, 1./2.) ;

            A_matrix[index_2d_to_1d(i, i_upwind, 
                                    system_rank, system_rank)] -= C_sat
                * sgn
                * dA 
                * dt
                / (V[i] * rho[i]);

            A_matrix[index_2d_to_1d(i+1, i_upwind, 
                                    system_rank, system_rank)] += C_sat
                * sgn
                * dA 
                * dt
                / (V[i+1] * rho[i+1]);
        }
    }

    for (int i=0 ; i<system_rank*system_rank; ++i)
    {
        if (!std::isfinite(A_matrix[i]))
        {
            printf("Non-finite value: A_matrix[i=%d]=%e\n", i, A_matrix[i]);
            assert(std::isfinite(A_matrix[i]));
        }
    }

    // Now it's time to call LAPACKE_dgbsv
    // here are the docs: http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html#gafa35ce1d7865b80563bbed6317050ad7

    lapack_int N = system_rank;
    lapack_int KL = 1;
    lapack_int KU = 1;
    lapack_int NRHS = 1;
    // lapack_int LDAB = (2*KL) + KU + 1
    lapack_int LDAB = N; // LAPACKE enforces LDAB >= N
    double *AB = new double[N*N]();
    lapack_int IPIV[N];
    lapack_int LDB = NRHS; // The documentation makes this sound like N, but it actually needs NRHS...

    double B[system_rank][1];
    for( int i=0; i<system_rank; ++i)
    {
        B[i][0] = e_int_old[i];
    }

    // it'll be less bug-prone if I just use their conventions,
    // and switch it back to 0 indexing at the end
    // And yes, I really do mean <= for the stopping conditions
    for( int j_fortran=1 ; j_fortran<=N ; ++j_fortran )
    {
        for( int i_fortran=std::max(1, j_fortran-KU) ; 
                 i_fortran<=std::min(N, j_fortran+KL) ;
                 ++i_fortran)
        {
            AB[index_2d_to_1d(KL + KU + 1 + i_fortran - j_fortran - 1, j_fortran-1,
                              N, N)] = A_matrix[index_2d_to_1d(i_fortran-1, j_fortran-1,
                                                               system_rank, system_rank)];
        }
    }


    if(verbose)
    {
        printf("about to call LAPACKE_dgbsv\n");
        fflush(stdout);
    }

    lapack_int INFO; // result status: 0 for success, <0 for (-INFO)th argument invalid; >0 and U(INFO, INFO)=0 for singular matrix and failed solution
    INFO = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,
        N, // the number of linear equations
        KL, // number of subdiagonals within the band
        KU, // number of superdiagonals within the band
        NRHS, // number of right hand sides (i.e. number of columns in B)
        AB, // on input, this is matrix A *in band storage*; on exit it's the factorization
        LDAB, // leading dimension of array AB (I think this is just N for me?)
        IPIV, // pivot indices; just an output I don't need (size N)
        *B, // B, on input it's B (NxNRHS); on exist if INFO=0, it's solution X
        LDB // The documentation makes this sound like N, but it actually needs NRHS...
        );

    if( INFO > 0 )
    {
        printf("INFO positive. (INFO = %d)\n", INFO);
        printf("That means it hit a singular matrix and failed.\n");
        fflush(stdout);
        assert(INFO==0);
    }
    else if( INFO < 0 )
    {
        printf("INFO negative. (INFO = %d)\n", INFO);
        printf("That means %dth input was invalid.\n", -INFO);
        fflush(stdout);
        assert(INFO==0);
    }


    for( int i=0 ; i<N ; ++i)
    {
        e_int_new[i] = B[i][0];
    }


    // ======== Verify post-conditions ========= //
    // #ifndef NDEBUG
    for( int i=0 ; i<num_cells ; ++i)
    {
        if( (i==0) && (verbose) )
        {
            printf("------ TEST in thermal_conduction_apply_linsolve() postconditions------- \n");
            printf("e_int_old[  i] = %e\n", e_int_old[i]);
            printf("e_int_guess[i] = %e\n", e_int_guess[i]);
            printf("e_int_new[  i] = %e\n", e_int_new[i]);

            printf("e_int_old[  i+1] = %e\n", e_int_old[  i+1]);
            printf("e_int_guess[i+1] = %e\n", e_int_guess[i+1]);
            printf("e_int_new[  i+1] = %e\n", e_int_new[  i+1]);

            const double mu_L = mu[i];
            const double mu_R = mu[i+1];

            const double mu_RR = mu[i+1];

            const double T_L_new = (mu_L * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_new[i];
            const double T_R_new = (mu_R * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_new[i+1];
            const double T_RR_new = (mu_RR * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_new[i+2];

            printf("T_L_new = %e\n", T_L_new);
            printf("T_R_new = %e\n", T_R_new);
            printf("T_RR_new = %e\n", T_RR_new);

            printf("dt = %e [s]\n", dt);
        }

        if( e_int_new[i] < 0 )
        {
            printf("------ ERROR in thermal_conduction_apply_linsolve() postconditions------- \n");
            printf("Negative specific internal energy within cell i=%d!\n", i);
            printf("e_int_old[  i] = %e\n", e_int_old[  i]);
            printf("e_int_guess[i] = %e\n", e_int_guess[i]);
            printf("e_int_new[  i] = %e\n", e_int_new[  i]);
            printf("dt = %e [s]\n", dt);

            fflush(stdout);

            assert(e_int_new[i] > 0);
        }
    }
    // #endif

    delete[] AB;
    delete[] A_matrix;
}

void spitzer_thermal_conduction_apply_linsolve( const double * e_int_old , 
                                const double * e_int_guess , 
                                double * e_int_new , 
                                const double * r ,
                                const double * V ,
                                const double * rho ,
                                const double * mu ,
                                const double dt ,
                                const int num_cells ,
                                const int verbose )
{

    // ============================================= //
    //
    //  Linearizes and solves the conduction equation for an implicit _guess_
    //  of e_int_guess. Returns a linearized solution e_int_new. 
    //  For Picard iteration, e_int_new becomes e_int_guess of the next iteration.
    //
    //  Inputs:
    //    - e_int_old, e_int_guess, e_int_new
    //        - array of internal energy per unit mass
    //        - has length num_cells
    //        - e_int_new is an output, and can be garbage on input
    //    - r  - array of cell-centered radii for each cell
    //    - V   - array of cell volumes
    //    - mu  - array of cell mean molecular weights
    //    - rho - array of cell densities
    //    - dt  - timestep
    //    - num_cells - the number of cells considered, not including guard cells
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - e_int_new completely overwritten with result
    //
    //  Notes:
    //    - The banded matrix solver we use is LAPACKE_dgbsv
    //      (This uses the c binding for the fortran LAPACK package)
    //    - For a good example of using LAPACK from c, and for the 
    //      indexing convention, check out:
    //      http://www.netlib.org/lapack/explore-html/d8/dd5/example___d_g_e_l_s__rowmajor_8c_source.html
    //    - I chose a banded diagonal solver when I was trying to solve a more
    //      complex problem. In principle, you could go back to using a simpler
    //      tri-diagonal solver, but I've got this debugged already
    //
    // ============================================= // 

    if(verbose)
    {
        printf("just entered thermal_conduction_apply_linsolve\n");
        fflush(stdout);
    }



    const int system_rank = num_cells;

    // Declare and initialize matrices so we can do:
    // A x = b
    // where b is `e_int_old`
    // and A is determined using `e_int_guess`

    double *A_matrix = new double[system_rank*system_rank];

    // initialize linear algebra containers
    for( int i=0 ; i<system_rank ; ++i )
    {
        for( int j=0 ; j<system_rank ; ++j )
        {
            A_matrix[index_2d_to_1d(i, j, system_rank, system_rank)] = 0.;
        }
    }

    if(verbose)
    {
        printf("adding identity matrix\n");
        fflush(stdout);
    }

    // add identity along diagonal
    for( int i=0 ; i<system_rank ; ++i )
    {
        A_matrix[index_2d_to_1d(i, i, system_rank, system_rank)] = 1.;
    }

    for( int i=0 ; i<(num_cells-1) ; ++i )
    {
        // this is interface i + 1/2
        const double r_interface = (r[i+1] + r[i]) / 2.;
        const double dA = get_dA(r_interface);
        const double dr = r[i+1] - r[i];

        const double mu_L = mu[i];
        const double mu_R = mu[i+1];

        const double e_int_L = e_int_guess[i];
        const double e_int_R = e_int_guess[i+1];

        const double T_L = (mu_L * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_L;
        const double T_R = (mu_R * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_R;

        // first, check if the temperatures are relatively close
        // then, check if there are effectively equal
        // if so, skip the conduction routine 
        //    (if you don't you might get silent divide-by-zero errors leading to NaN energy flux)
        if (std::abs(std::log10(T_L/T_R)) < 1)
        {
            if (std::abs(1 - (T_L/T_R)) < .01)
            {
                // effectively no temperature gradient; don't apply conduction
                // otherwise you get silent divide-by-zero errors, leading to 
                // fluxes that are NaN's
                continue;
            }
        }

        // okay, so now we need to know if it's going to be saturated or unsaturated.
        // The only real way to do this is actually calculate the energy flux.

        double rho_upwind;
        double e_int_upwind;

        int sgn; // sign of dT/dr deriv; (+1 if energy is transfered left, -1 if energy transfered right)
        int i_upwind;
        int i_downwind;
        if( e_int_L > e_int_R )
        {   
            i_upwind   = i;
            i_downwind = i+1;
            sgn        = -1;
            // left cell is hotter and therefore "upwind"
            rho_upwind = rho[i_upwind];
            e_int_upwind = e_int_guess[i_upwind];
        }
        else
        {
            i_upwind   = i + 1;
            i_downwind = i;
            sgn        = +1;
            // right cell is hotter and therefore "upwind"
            rho_upwind = rho[i_upwind];
            e_int_upwind = e_int_guess[i_upwind];
        }

        const double mu_average = (mu_L + mu_R) / 2.;
        const double e_int_average = (e_int_L + e_int_R) / 2.;


        const double electron_charge = 4.8032047e-10; // Fr
        const double electron_mass = 9.109383e-28;

        const double C_spitzer = std::pow(e_int_average, 5./2.) 
                * std::pow( mu_average * m_proton * (GAMMA_LAW-1) , 7./2.)
                / (std::pow(electron_mass, .5) * std::pow(electron_charge, 4.)) ;





        A_matrix[index_2d_to_1d(i, i, 
                                system_rank, system_rank)] += C_spitzer
            * dA
            * dt
            / (V[i] * rho[i])
            / dr;

        A_matrix[index_2d_to_1d(i, i+1, 
                                system_rank, system_rank)] -= C_spitzer
            * dA
            * dt
            / (V[i] * rho[i])
            / dr;


        A_matrix[index_2d_to_1d(i+1, i, 
                                system_rank, system_rank)] -= C_spitzer
            * dA
            * dt
            / (V[i+1] * rho[i+1])
            / dr;


        A_matrix[index_2d_to_1d(i+1, i+1, 
                                system_rank, system_rank)] += C_spitzer
            * dA
            * dt
            / (V[i+1] * rho[i+1])
            / dr;
    }


    for (int i=0 ; i<system_rank*system_rank; ++i)
    {
        if (!std::isfinite(A_matrix[i]))
        {
            printf("Non-finite value: A_matrix[i=%d]=%e\n", i, A_matrix[i]);
            assert(std::isfinite(A_matrix[i]));
        }
    }

    // Now it's time to call LAPACKE_dgbsv
    // here are the docs: http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html#gafa35ce1d7865b80563bbed6317050ad7

    lapack_int N = system_rank;
    lapack_int KL = 1;
    lapack_int KU = 1;
    lapack_int NRHS = 1;
    // lapack_int LDAB = (2*KL) + KU + 1
    lapack_int LDAB = N; // LAPACKE enforces LDAB >= N
    double *AB = new double[N*N]();
    lapack_int IPIV[N];
    lapack_int LDB = NRHS; // The documentation makes this sound like N, but it actually needs NRHS...

    double B[system_rank][1];
    for( int i=0; i<system_rank; ++i)
    {
        B[i][0] = e_int_old[i];
    }

    // it'll be less bug-prone if I just use their conventions,
    // and switch it back to 0 indexing at the end
    // And yes, I really do mean <= for the stopping conditions
    for( int j_fortran=1 ; j_fortran<=N ; ++j_fortran )
    {
        for( int i_fortran=std::max(1, j_fortran-KU) ; 
                 i_fortran<=std::min(N, j_fortran+KL) ;
                 ++i_fortran)
        {
            AB[index_2d_to_1d(KL + KU + 1 + i_fortran - j_fortran - 1, j_fortran-1,
                              N, N)] = A_matrix[index_2d_to_1d(i_fortran-1, j_fortran-1,
                                                               system_rank, system_rank)];
        }
    }


    if(verbose)
    {
        printf("about to call LAPACKE_dgbsv\n");
        fflush(stdout);
    }

    lapack_int INFO; // result status: 0 for success, <0 for (-INFO)th argument invalid; >0 and U(INFO, INFO)=0 for singular matrix and failed solution
    INFO = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,
        N, // the number of linear equations
        KL, // number of subdiagonals within the band
        KU, // number of superdiagonals within the band
        NRHS, // number of right hand sides (i.e. number of columns in B)
        AB, // on input, this is matrix A *in band storage*; on exit it's the factorization
        LDAB, // leading dimension of array AB (I think this is just N for me?)
        IPIV, // pivot indices; just an output I don't need (size N)
        *B, // B, on input it's B (NxNRHS); on exist if INFO=0, it's solution X
        LDB // The documentation makes this sound like N, but it actually needs NRHS...
        );

    if( INFO > 0 )
    {
        printf("INFO positive. (INFO = %d)\n", INFO);
        printf("That means it hit a singular matrix and failed.\n");
        fflush(stdout);
        assert(INFO==0);
    }
    else if( INFO < 0 )
    {
        printf("INFO negative. (INFO = %d)\n", INFO);
        printf("That means %dth input was invalid.\n", -INFO);
        fflush(stdout);
        assert(INFO==0);
    }


    for( int i=0 ; i<N ; ++i)
    {
        e_int_new[i] = B[i][0];
    }


    // ======== Verify post-conditions ========= //
    // #ifndef NDEBUG
    for( int i=0 ; i<num_cells ; ++i)
    {
        if( (i==0) && (verbose) )
        {
            printf("------ TEST in thermal_conduction_apply_linsolve() postconditions------- \n");
            printf("e_int_old[  i] = %e\n", e_int_old[i]);
            printf("e_int_guess[i] = %e\n", e_int_guess[i]);
            printf("e_int_new[  i] = %e\n", e_int_new[i]);

            printf("e_int_old[  i+1] = %e\n", e_int_old[  i+1]);
            printf("e_int_guess[i+1] = %e\n", e_int_guess[i+1]);
            printf("e_int_new[  i+1] = %e\n", e_int_new[  i+1]);

            const double mu_L = mu[i];
            const double mu_R = mu[i+1];

            const double mu_RR = mu[i+1];

            const double T_L_new = (mu_L * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_new[i];
            const double T_R_new = (mu_R * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_new[i+1];
            const double T_RR_new = (mu_RR * m_proton * (GAMMA_LAW-1) / k_boltzmann) * e_int_new[i+2];

            printf("T_L_new = %e\n", T_L_new);
            printf("T_R_new = %e\n", T_R_new);
            printf("T_RR_new = %e\n", T_RR_new);

            printf("dt = %e [s]\n", dt);
        }

        if( e_int_new[i] < 0 )
        {
            printf("------ ERROR in thermal_conduction_apply_linsolve() postconditions------- \n");
            printf("Negative specific internal energy within cell i=%d!\n", i);
            printf("e_int_old[  i] = %e\n", e_int_old[  i]);
            printf("e_int_guess[i] = %e\n", e_int_guess[i]);
            printf("e_int_new[  i] = %e\n", e_int_new[  i]);
            printf("dt = %e [s]\n", dt);

            fflush(stdout);

            assert(e_int_new[i] > 0);
        }
    }
    // #endif

    delete[] AB;
    delete[] A_matrix;
}

void turbulent_diffusion_apply_linsolve( const double * Y_old ,
                                         const double * Y_guess ,
                                         double * Y_new ,
                                         const double * r ,
                                         const double * Vol ,
                                         const double * v ,
                                         const double dt ,
                                         const double C_turbulent_diffusion ,
                                         const int num_cells ,
                                         const int verbose,
                                         const bool allow_negatives )
{

    // ============================================= //
    //
    //  Note: this applies diffusion on a conservative variable _per unit volume_.
    //  e.g. it applies to metal density, not metallicity.
    //
    //  Also note, this doesn't actually need to be iterated on! Since the 
    //  matrix _never changes_ (since I don't update v), it only needs one
    //  linear solve.
    //
    //  Linearizes and solves the turbulent diffusion equation for an implicit 
    //  _guess_ of Y_guess. Note: this part doesn't care about what `Y` is,
    //  as long as it's a conservative variable per unit volume.
    //  But in practice, we'll only be solving the mass flux for
    //  turbulent diffusion, and then assume everything else traces the mass.
    //
    //  Returns a linearized solution Y_new. 
    //  For Picard iteration, Y_new becomes Y_guess of the next iteration.
    // 
    //  Within my notes I use `A` as the conservative variable, but that'll get
    //  confusing because I write the linearized matrix as `A` below. So `Y`
    //  seemed to have the least amount of conflicts.
    //
    //  Inputs:
    //    - Y_old, Y_guess, Y_new
    //        - array of conservative per unit *volume*
    //        - has length num_cells
    //        - Y_new is an output, and can be garbage on input
    //    - r   - array of cell-centered radii for each cell
    //    - Vol - array of cell volumes
    //    - v - array of cell fluid velocities
    //    - dt  - timestep
    //    - C_turbulent_diffusion - scales the location-variable diffusion coeff.
    //                              ( D = C |S| h^2 )
    //    - num_cells - the number of cells considered, not including guard cells
    //    - verbose
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - Y_new completely overwritten with result
    //
    //  Notes:
    //    - The banded matrix solver we use is LAPACKE_dgbsv
    //      (This uses the c binding for the fortran LAPACK package)
    //    - For a good example of using LAPACK from c, and for the 
    //      indexing convention, check out:
    //      http://www.netlib.org/lapack/explore-html/d8/dd5/example___d_g_e_l_s__rowmajor_8c_source.html
    //    - I chose a banded diagonal solver when I was trying to solve a more
    //      complex problem. In principle, you could go back to using a simpler
    //      tri-diagonal solver, but I've got this debugged already
    //
    // ============================================= // 


    if(verbose)
    {
        printf("just entered turbulent_diffusion_apply_linsolve\n");
        fflush(stdout);
    }


    const int system_rank = num_cells;

    // Declare and initialize matrices so we can do:
    // A x = b
    // where b is `Y_old`
    // and A is determined using `Y_guess`

    double *A_matrix = new double[system_rank*system_rank];

    // initialize linear algebra containers
    for( int i=0 ; i<system_rank ; ++i )
    {
        for( int j=0 ; j<system_rank ; ++j )
        {
            A_matrix[index_2d_to_1d(i, j, system_rank, system_rank)] = 0.;
        }
    }

    if(verbose)
    {
        printf("adding identity matrix\n");
        fflush(stdout);
    }

    // add identity along diagonal
    for( int i=0 ; i<system_rank ; ++i )
    {
        A_matrix[index_2d_to_1d(i, i, system_rank, system_rank)] = 1.;
    }



    for( int i=0 ; i<(num_cells-1) ; ++i )
    {

        const double b = turbulent_diffusion_calculate_b(r[i], r[i+1],
                                                         v[i], v[i+1], 
                                                         C_turbulent_diffusion);

        A_matrix[index_2d_to_1d(i, i, 
                                system_rank, system_rank)] += b
            * dt
            / Vol[i];

        A_matrix[index_2d_to_1d(i, i+1, 
                                system_rank, system_rank)] -= b
            * dt
            / Vol[i];


        A_matrix[index_2d_to_1d(i+1, i, 
                                system_rank, system_rank)] -= b
            * dt
            / Vol[i+1];


        A_matrix[index_2d_to_1d(i+1, i+1, 
                                system_rank, system_rank)] += b
            * dt
            / Vol[i+1];

    }

    for (int i=0 ; i<system_rank*system_rank; ++i)
    {
        if (!std::isfinite(A_matrix[i]))
        {
            printf("Non-finite value: A_matrix[i=%d]=%e\n", i, A_matrix[i]);
            assert(std::isfinite(A_matrix[i]));
        }
    }

    // Now it's time to call LAPACKE_dgbsv
    // here are the docs: http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html#gafa35ce1d7865b80563bbed6317050ad7

    lapack_int N = system_rank;
    lapack_int KL = 1;
    lapack_int KU = 1;
    lapack_int NRHS = 1;
    // lapack_int LDAB = (2*KL) + KU + 1
    lapack_int LDAB = N; // LAPACKE enforces LDAB >= N
    double *AB = new double[N*N]();
    lapack_int IPIV[N];
    lapack_int LDB = NRHS; // The documentation makes this sound like N, but it actually needs NRHS...

    double B[system_rank][1];
    for( int i=0; i<system_rank; ++i)
    {
        B[i][0] = Y_old[i];
    }

    // it'll be less bug-prone if I just use their conventions,
    // and switch it back to 0 indexing at the end
    // And yes, I really do mean <= for the stopping conditions
    for( int j_fortran=1 ; j_fortran<=N ; ++j_fortran )
    {
        for( int i_fortran=std::max(1, j_fortran-KU) ; 
                 i_fortran<=std::min(N, j_fortran+KL) ;
                 ++i_fortran)
        {
            AB[index_2d_to_1d(KL + KU + 1 + i_fortran - j_fortran - 1, j_fortran-1,
                              N, N)] = A_matrix[index_2d_to_1d(i_fortran-1, j_fortran-1,
                                                               system_rank, system_rank)];
        }
    }

    if(verbose)
    {
        printf("about to call LAPACKE_dgbsv\n");
        fflush(stdout);
    }

    lapack_int INFO; // result status: 0 for success, <0 for (-INFO)th argument invalid; >0 and U(INFO, INFO)=0 for singular matrix and failed solution
    INFO = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,
        N, // the number of linear equations
        KL, // number of subdiagonals within the band
        KU, // number of superdiagonals within the band
        NRHS, // number of right hand sides (i.e. number of columns in B)
        AB, // on input, this is matrix A *in band storage*; on exit it's the factorization
        LDAB, // leading dimension of array AB (I think this is just N for me?)
        IPIV, // pivot indices; just an output I don't need (size N)
        *B, // B, on input it's B (NxNRHS); on exist if INFO=0, it's solution X
        LDB // The documentation makes this sound like N, but it actually needs NRHS...
        );

    if( INFO > 0 )
    {
        printf("INFO positive. (INFO = %d)\n", INFO);
        printf("That means it hit a singular matrix and failed.\n");
        fflush(stdout);
        assert(INFO==0);
    }
    else if( INFO < 0 )
    {
        printf("INFO negative. (INFO = %d)\n", INFO);
        printf("That means %dth input was invalid.\n", -INFO);
        fflush(stdout);
        assert(INFO==0);
    }


    for( int i=0 ; i<N ; ++i)
    {
        Y_new[i] = B[i][0];
    }


    // ======== Verify post-conditions ========= //
    // #ifndef NDEBUG
    for( int i=0 ; i<num_cells ; ++i)
    {
        if( (i==0) && (verbose) )
        {
            printf("------ TEST in turbulent_diffusion_apply_linsolve() postconditions------- \n");
            printf("Y_old[  i] = %e\n", Y_old[i]);
            printf("Y_guess[i] = %e\n", Y_guess[i]);
            printf("Y_new[  i] = %e\n", Y_new[i]);

            printf("Y_old[  i+1] = %e\n", Y_old[  i+1]);
            printf("Y_guess[i+1] = %e\n", Y_guess[i+1]);
            printf("Y_new[  i+1] = %e\n", Y_new[  i+1]);

            printf("dt = %e [s]\n", dt);
        }

        if( !allow_negatives )
        {
            if( Y_new[i] < 0 )
            {
                printf("------ ERROR in turbulent_diffusion_apply_linsolve() postconditions------- \n");
                printf("Negative Y within cell i=%d!\n", i);
                printf("Y_old[  i] = %e\n", Y_old[  i]);
                printf("Y_guess[i] = %e\n", Y_guess[i]);
                printf("Y_new[  i] = %e\n", Y_new[  i]);
                printf("dt = %e [s]\n", dt);

                fflush(stdout);

                assert(Y_new[i] > 0);
            }
        }

    }
    // #endif

    delete[] AB;
    delete[] A_matrix;
}

double turbulent_diffusion_calculate_b(
    const double r_left ,
    const double r_right ,
    const double v_left ,
    const double v_right ,
    const double C_turbulent_diffusion
     )
{
    // ============================================= //
    //
    //  Calculates the prefactor for artificial turbulent diffusion.
    //
    //  This prefactor is not a commonly recognized definition; it's just a way
    //  to make my math a little simpler, and allows the core turbulent diffusion
    //  prescription to be reused in multiple areas without having to copy and
    //  and paste.
    // 
    //  This prefactor is defined such that the mass fluxes in spherical 
    //  symmetry are simply: 
    //      mass flux_{i + 1/2} = - b * (density_{i+1} - density_{i})
    //  The general ideal is described at:
    //    https://github.com/egentry/clustered_SNe/wiki/Turbulent-diffusion-prescription
    //  I haven't yet pushed my notebook going more into the math, but it
    //    *might* be visible here:
    //    https://www.dropbox.com/s/xowqsculq4t6yzf/turbulent%20diffusion%20write%20up.ipynb?dl=0
    //
    //
    //  Inputs:
    //    - r_left, r_right
    //        - cell-centered radii for the left and right cells ("right" means
    //          larger radius)
    //    - v_left, v_right
    //        - fluid velocities for the left and right cells
    //    - C_turbulent_diffusion - scales the location-variable diffusion coeff.
    //                              ( D = C |S| h^2 )
    //
    //  Returns:
    //    - turbulent_diffusion_b
    //      - should always be non-negative
    //      - may be zeroed out for small values to avoid numerical noise
    //
    // ============================================= // 

        // this is interface i + 1/2
        const double r_interface = (r_right + r_left) / 2.;
        const double dr = r_right - r_left;

        // // // Create diffusion coefficient, D = C |S| h^2

        const double dv_dr = (v_right - v_left) / dr;
        const double v_interface = (v_right + v_left) / 2;

        if(std::abs(dv_dr - (v_interface / r_interface)) < 1. / pc_unit )
        {
            // | dv/dr  - (v/r)| < 1 cm/s / pc
            // velocity term negligible; don't waste time on noise
            return 0;
        }

        const double S_mag = std::sqrt(2./3.) 
                             * std::abs(dv_dr - (v_interface/r_interface));
        const double D = C_turbulent_diffusion * S_mag * std::pow(dr, 2.);
        const double turbulent_diffusion_b = 4. * M_PI * std::pow(r_interface, 2.) * D / dr;

        // verify postconditions
        #ifndef NDEBUG
        if( turbulent_diffusion_b < 0 )
        {
            printf("Value Error: turbulent_diffusion_b = %e but should be non-negative\n" ,
                turbulent_diffusion_b);
            fflush(stdout);
            assert(turbulent_diffusion_b >= 0);
        }

        if( !std::isfinite(turbulent_diffusion_b) )
        {
            printf("Value Error: turbulent_diffusion_b =  %e for interface \n", 
                turbulent_diffusion_b);
            fflush(stdout);
            assert(std::isfinite(turbulent_diffusion_b));
        }
        #endif

        return turbulent_diffusion_b;
}


void turbulent_diffusion_get_fluxes( const double * densities , 
                                           double * mass_fluxes ,
                                     const double * r ,
                                     const double * v ,
                                     const double C_turbulent_diffusion ,
                                     const int num_cells ,
                                     const int verbose )
{

    // ============================================= //
    //
    //  Calculates the mass fluxes between all cells using the artificial
    //  turbulent diffusion prescription given an array of densities.
    // 
    //  This will look a little different from `turbulent_diffusion_apply_linsolve`
    //  because this will look more like the explicit approach to dealing with
    //  diffusion. The key difference from a true explicit approach is that
    //  here we take a field of densities which _are not_ the current densities.
    //  Instead, we should be using the densities _which will result from the
    //  calculated fluxes_ (making this approach implicit).
    //
    //  The reason why we aren't simply using the already-calculated new densities
    //  is that we to make sure those mass fluxes also carry corresponding
    //  fluxes of the other conservatives, proportional to the mass flux. 
    //  I.e. the energy flux should be:
    //     energy_flux = cell_energy * (mass_flux / cell_mass)
    //
    //  Inputs:
    //    - densities
    //        - array of densities, not including ghost zones
    //        - has length num_cells
    //    - mass_fluxes
    //        - array of _output_ mass fluxes
    //        - has length of num_cells-1
    //        - can be garbage on input
    //        - if mass_fluxes[i] is positive, means flux from cell i to i+1
    //          if negative, means flux from cell i+1 to i
    //    - r   - array of cell-centered radii for each cell
    //    - v - array of cell fluid velocities
    //    - C_turbulent_diffusion - scales the location-variable diffusion coeff.
    //                              ( D = C |S| h^2 )
    //    - num_cells - the number of cells considered, not including guard cells
    //    - verbose
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - mass_fluxes completely overwritten with result
    //
    // ============================================= // 


    if(verbose)
    {
        printf("just entered turbulent_diffusion_get_fluxes\n");
        fflush(stdout);
    }
 
    for( int i=0 ; i<(num_cells-1) ; ++i )
    {
        const double b = turbulent_diffusion_calculate_b(r[i], r[i+1],
                                                         v[i], v[i+1],
                                                         C_turbulent_diffusion );
        mass_fluxes[i] = - b * (densities[i+1] - densities[i]);
    }
}

void U_i_to_prim( const double * U , double * prim , const int i )
{

    const double rho = U[(4*i) + RHO];
    const double v   = U[(4*i) + VRR] / rho;
    const double P = (GAMMA_LAW-1)*(U[(4*i)+PPP] - .5*std::pow(U[(4*i)+VRR],2.) / U[(4*i)+RHO]);
    const double Z   = U[(4*i) + ZZZ] / rho;

    prim[RHO] = rho;
    prim[VRR] = v;
    prim[PPP] = P;
    prim[ZZZ] = Z;
}

double max_norm( const double * A , const double * B, const int size )
{

    // ============================================= //
    //
    //  Calculate the L-infinity (max) norm between two arrays
    //
    //  Inputs:
    //    - A, B - vectors of equal size
    //    - size - the size of the arrays
    //
    //  Returns:
    //    - max_distance - the max difference between corresponding elements
    //
    //  Side effects:
    //    None
    //
    // ============================================= // 

    double max_distance = 0;
    for( int i=0; i<size ; ++i )
    {
        const double distance = std::abs(A[i] - B[i]);
        if( distance > max_distance )
        {
            max_distance = distance;
        }
    }

    return max_distance;
}

void thermal_conduction_implicit( struct domain * theDomain , const double dt )
{

    // ============================================= //
    //
    //  Calculate and add heat fluxes from *physical* conduction.
    //
    //  Inputs:
    //    - theDomain    - the standard domain struct used throughout
    //    - dt           - timestep [s]
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - affects cons for each cell.
    //
    //  Notes:
    //    - The banded matrix solver we use is LAPACKE_dgbsv
    //      (This uses the c binding for the fortran LAPACK package)
    //    - For a good example of using LAPACK from c, and for the 
    //      indexing convention, check out:
    //      http://www.netlib.org/lapack/explore-html/d8/dd5/example___d_g_e_l_s__rowmajor_8c_source.html
    //
    // ============================================= // 

    // printf("Entered thermal_conduction_implicit\n");
    // fflush(stdout);


    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;

    const int num_cells = Nr-Ng-Ng;

    double e_int_old[num_cells];
    double e_int_guess[num_cells];
    double e_int_new[num_cells];

    double r[num_cells];
    double V[num_cells];
    double rho[num_cells];
    double mu[num_cells];

    for( int i_old_arrays=Ng ; i_old_arrays<Nr-Ng ; ++i_old_arrays )
    {
        const struct cell * c = &(theCells[i_old_arrays]);

        const int i_new_arrays = i_old_arrays - 1;

        r[i_new_arrays] = c->riph - (c->dr/2.); // cell *centered* radius
        V[i_new_arrays] = get_dV(c->riph, c->riph - c->dr); // cell volume

        rho[i_new_arrays] = c->prim[RHO];
        mu[i_new_arrays] = get_mean_molecular_weight(c->prim[ZZZ]);

        e_int_old[i_new_arrays] = c->prim[PPP] / c->prim[RHO] / (GAMMA_LAW-1);
    }

    // initialize my first e_int_guess using the previous conditions
    for( int i=0 ; i<num_cells ; ++i )
    {
        e_int_guess[i] = e_int_old[i];
    }

    std::function<void(const double *, 
                       const double *, 
                             double *)> lambda_apply_linsolve = [&] (const double * values_old,
                                            const double * values_guess,
                                            double * values_new) {
        thermal_conduction_apply_linsolve(values_old, values_guess, values_new, 
                                  r, V, rho, mu, dt, num_cells, 0 );
        // spitzer_thermal_conduction_apply_linsolve(values_old, values_guess, values_new, 
                                  // r, V, rho, mu, dt, num_cells, 0 );

    };

    if(theDomain->theParList.with_anderson_acceleration)
    {
        Anderson_acceleration(e_int_old, e_int_guess, e_int_new,
                              num_cells, lambda_apply_linsolve, 
                              "thermal_conduction_implicit");
    }
    else
    {
        Picard_iteration(e_int_old, e_int_guess, e_int_new,
                         num_cells, lambda_apply_linsolve, 
                         "thermal_conduction_implicit");
    }


    // Update cons using the difference in e_int
    double dE_total = 0;
    for( int i_linsolve=0 ; i_linsolve<num_cells ; ++i_linsolve )
    {
        const int i_cells = i_linsolve + 1; 
        struct cell * c = &(theCells[i_cells]);

        const double mass = rho[i_linsolve] * V[i_linsolve];
        const double de_int = e_int_new[i_linsolve] - e_int_old[i_linsolve];
        const double dE = mass * de_int;

        c->cons[TAU] += dE;
        dE_total += dE;
    }
    theDomain->energy_added_by_thermal_conduction += dE_total;


    // ======== Verify post-conditions ========= //
    // #ifndef NDEBUG
    if(std::abs(dE_total) > (1e51 * 1e-10))
    {
        double E_total = 0;
        for( int i_linsolve=0 ; i_linsolve<num_cells ; ++i_linsolve )
        {
            const int i_cells = i_linsolve + 1; 
            const struct cell * c = &(theCells[i_cells]);
            E_total += c->cons[TAU];
        }
        printf("------ ERROR in thermal_conduction_implicit() postconditions------- \n");
        printf("dE_total outside of threshold! \n");
        printf("dE_total  = %e \n", dE_total);
        printf("E_total   = %e \n", E_total);

        assert(std::abs(dE_total) < (1e51 * 1e-10));   
    }

    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        struct cell * c = &(theCells[i]);
        if( c->cons[DDD] <= 0 )
        {
            printf("------ ERROR in thermal_conduction_implicit() postconditions------- \n");
            printf("Cell i=%d has negative mass! \n", i);
            printf("c->cons[DDD] = %e \n", c->cons[DDD]);
            printf("c->cons[SRR] = %e \n", c->cons[SRR]);
            printf("c->cons[TAU] = %e \n", c->cons[TAU]);
            printf("c->cons[ZZZ] = %e \n", c->cons[ZZZ]);

            assert(c->cons[DDD] > 0);
        }

        if( c->cons[TAU] <= 0 )
        {
            printf("------ ERROR in thermal_conduction_implicit() postconditions------- \n");
            printf("Cell i=%d has negative energy! \n", i);
            printf("c->cons[DDD] = %e \n", c->cons[DDD]);
            printf("c->cons[SRR] = %e \n", c->cons[SRR]);
            printf("c->cons[TAU] = %e \n", c->cons[TAU]);
            printf("c->cons[ZZZ] = %e \n", c->cons[ZZZ]);

            assert(c->cons[TAU] > 0);
        }

        if( c->cons[ZZZ] <= 0 )
        {
            printf("------ ERROR in thermal_conduction_implicit() postconditions------- \n");
            printf("Cell i=%d has negative metal content! \n", i);
            printf("c->cons[DDD] = %e \n", c->cons[DDD]);
            printf("c->cons[SRR] = %e \n", c->cons[SRR]);
            printf("c->cons[TAU] = %e \n", c->cons[TAU]);
            printf("c->cons[ZZZ] = %e \n", c->cons[ZZZ]);

            assert(c->cons[ZZZ] > 0);
        }
    }
    // #endif

}




void turbulent_diffusion_implicit( struct domain * theDomain , 
                                   const double dt )
{
    // ============================================= //
    //
    //  Calculate and apply mass-tracing fluxes due to unresolved
    //  turbulent mixing.
    //
    //  This'll first calculate the mass fluxes implicitly, then will scale those
    //  fluxes to the other conserved variables.
    //
    //  The general ideal is described at:
    //    https://github.com/egentry/clustered_SNe/wiki/Turbulent-diffusion-prescription
    //  I haven't yet pushed my notebook going more into the math, but it
    //    *might* be visible here:
    //    https://www.dropbox.com/s/xowqsculq4t6yzf/turbulent%20diffusion%20write%20up.ipynb?dl=0
    //
    //  Inputs:
    //    - theDomain    - the standard domain struct used throughout
    //    - dt           - timestep [s]
    //
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - affects cons for each cell.
    //
    //  Notes:
    //    - The banded matrix solver we use is LAPACKE_dgbsv
    //      (This uses the c binding for the fortran LAPACK package)
    //    - For a good example of using LAPACK from c, and for the 
    //      indexing convention, check out:
    //      http://www.netlib.org/lapack/explore-html/d8/dd5/example___d_g_e_l_s__rowmajor_8c_source.html
    //
    // ============================================= // 

    // printf("Entered turbulent_diffusion_implicit\n");
    // fflush(stdout);


    struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;

    const int num_cells = Nr-Ng-Ng;

    double mass_density_old[num_cells];
    double mass_density_guess[num_cells];
    double mass_density_new[num_cells];

    double r[num_cells];
    double Vol[num_cells];
    double v[num_cells];

    for( int i_old_arrays=Ng ; i_old_arrays<Nr-Ng ; ++i_old_arrays )
    {
        const struct cell * c = &(theCells[i_old_arrays]);

        const int i_new_arrays = i_old_arrays - 1;

        r[i_new_arrays] = c->riph - (c->dr/2.); // cell *centered* radius
        Vol[i_new_arrays] = get_dV(c->riph, c->riph - c->dr); // cell volume

        v[i_new_arrays] = c->prim[VRR];

        mass_density_old[i_new_arrays] = c->prim[RHO];
    }

    // initialize my first mass_density_guess using the previous conditions
    for( int i=0 ; i<num_cells ; ++i )
    {
        mass_density_guess[i] = mass_density_old[i];
    }

    std::function<void(const double *, 
                       const double *, 
                             double *)> lambda_apply_linsolve = [&] (const double * values_old,
                                            const double * values_guess,
                                            double * values_new) {
        turbulent_diffusion_apply_linsolve(values_old, values_guess, values_new,
                                           r, Vol, v, dt, 
                                           theDomain->theParList.C_turbulent_diffusion,
                                           num_cells, 0 );
    };

    if(theDomain->theParList.with_anderson_acceleration)
    {
        Anderson_acceleration(mass_density_old, mass_density_guess, 
                              mass_density_new,
                              num_cells, lambda_apply_linsolve,
                              "turbulent_diffusion_implicit");
    }
    else
    {
        Picard_iteration(mass_density_old, mass_density_guess, mass_density_new,
                         num_cells, lambda_apply_linsolve,
                         "turbulent_diffusion_implicit");
    }

    double mass_fluxes[num_cells-1];
    turbulent_diffusion_get_fluxes(mass_density_new, mass_fluxes,
                                   r, v, 
                                   theDomain->theParList.C_turbulent_diffusion,
                                   num_cells);


    // Update cons using mass-tracing fluxes
    for( int i_linsolve=0 ; i_linsolve<num_cells-1 ; ++i_linsolve )
    {
        const int i_cells = i_linsolve + 1; 
        struct cell * c_left = &(theCells[i_cells]);
        struct cell * c_right = &(theCells[i_cells+1]);

        struct cell * c_upwind;
        struct cell * c_downwind;

        double d_mass; // will be non-negative / unsigned!

        if(mass_fluxes[i_linsolve] > 0)
        {
            c_upwind   = c_left;
            c_downwind = c_right;
            d_mass = mass_fluxes[i_linsolve] * dt;
        }   
        else
        {
            c_upwind   = c_right;
            c_downwind = c_left;
            d_mass = - mass_fluxes[i_linsolve] * dt;
        }

        const double d_momentum = d_mass *(c_upwind->cons[SRR] / c_upwind->cons[DDD]);
        const double d_energy   = d_mass *(c_upwind->cons[TAU] / c_upwind->cons[DDD]);
        const double d_metals   = d_mass *(c_upwind->cons[ZZZ] / c_upwind->cons[DDD]);

        c_upwind->cons[DDD]   -= d_mass;
        c_downwind->cons[DDD] += d_mass;

        c_upwind->cons[SRR]   -= d_momentum;
        c_downwind->cons[SRR] += d_momentum;

        c_upwind->cons[TAU]   -= d_energy;
        c_downwind->cons[TAU] += d_energy;

        c_upwind->cons[ZZZ]   -= d_metals;
        c_downwind->cons[ZZZ] += d_metals;
    }


    // ======== Verify post-conditions ========= //
    // #ifndef NDEBUG
    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        struct cell * c = &(theCells[i]);
        if( c->cons[DDD] <= 0 )
        {
            printf("------ ERROR in turbulent_diffusion_implicit() postconditions------- \n");
            printf("Cell i=%d has negative mass! \n", i);
            printf("c->cons[DDD] = %e \n", c->cons[DDD]);
            printf("c->cons[SRR] = %e \n", c->cons[SRR]);
            printf("c->cons[TAU] = %e \n", c->cons[TAU]);
            printf("c->cons[ZZZ] = %e \n", c->cons[ZZZ]);

            assert(c->cons[DDD] > 0);
        }

        if( c->cons[TAU] <= 0 )
        {
            printf("------ ERROR in turbulent_diffusion_implicit() postconditions------- \n");
            printf("Cell i=%d has negative energy! \n", i);
            printf("c->cons[DDD] = %e \n", c->cons[DDD]);
            printf("c->cons[SRR] = %e \n", c->cons[SRR]);
            printf("c->cons[TAU] = %e \n", c->cons[TAU]);
            printf("c->cons[ZZZ] = %e \n", c->cons[ZZZ]);

            assert(c->cons[TAU] > 0);
        }

        if( c->cons[ZZZ] <= 0 )
        {
            printf("------ ERROR in turbulent_diffusion_implicit() postconditions------- \n");
            printf("Cell i=%d has negative metal content! \n", i);
            printf("c->cons[DDD] = %e \n", c->cons[DDD]);
            printf("c->cons[SRR] = %e \n", c->cons[SRR]);
            printf("c->cons[TAU] = %e \n", c->cons[TAU]);
            printf("c->cons[ZZZ] = %e \n", c->cons[ZZZ]);

            assert(c->cons[ZZZ] > 0);
        }
    }
    // #endif

}

void subgrid_thermal_conduction( struct cell * c , const double dt )
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
    if (c->multiphase)
    {
        verify_multiphase_conditions(c, "subgrid_thermal_conduction()", "pre");

        // if (E_int_from_cons(c->cons_hot) < 0)
        // {
        //     printf("------ERROR in subgrid_thermal_conduction() preconditions------- \n");
        //     printf("hot gas *internal* energy less than 0\n");
        //     printf("E_int (hot)  = %e \n", E_int_from_cons(c->cons_hot));
        //     printf("c->cons_hot[TAU]   = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons[TAU]       = %e \n", c->cons[TAU]);

        //     assert( E_int_from_cons(c->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(c->cons_cold) < 0)
        // {
        //     printf("------ERROR in subgrid_thermal_conduction() preconditions------- \n");
        //     printf("cold gas *internal* energy less than 0\n");
        //     printf("E_int (cold)  = %e \n", E_int_from_cons(c->cons_cold));
        //     printf("c->cons_cold[TAU]   = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]        = %e \n", c->cons[TAU]);

        //     assert( E_int_from_cons(c->cons_cold) > 0 );
        // }

        // if ( (c->cons_hot[TAU] / c->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in subgrid_thermal_conduction() preconditions------- \n");
        //     printf("hot gas energy greater than total energy\n");
        //     printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

        //     assert( (c->cons_hot[TAU] / c->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( (c->cons_cold[TAU] / c->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in subgrid_thermal_conduction() preconditions------- \n");
        //     printf("cold gas energy greater than total energy\n");
        //     printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]      = %e \n", c->cons[TAU]);

        //     assert( (c->cons_cold[TAU] / c->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( std::abs(1-( (c->cons_cold[TAU] + c->cons_hot[TAU])/c->cons[TAU])) > rel_tol)
        // {
        //     printf("------ERROR in subgrid_thermal_conduction() preconditions------- \n");
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
    // // get the effective radius of the cold phase, assuming it's one sphere
    const double s_cold = std::pow((3. * c->V_cold / (4. * M_PI)), 1./3.);
    // // get the surface area of a sphere of that size
    const double A_cold = 4 * M_PI * std::pow(s_cold, 2);

    const double phi_s = 1.1; // see note near Eq. 8 of Cowie & McKee 1977
    const double c_s_cold = std::sqrt(std::abs(GAMMA_LAW*c->prim_cold[PPP]/c->prim_cold[RHO]));

    const double q = 5 * phi_s * c_s_cold * c->prim_cold[PPP];
    const double dM_b_dt_sat = q * (A_cold) / (c_s_cold*c_s_cold); // it's the colder cell's energy which is transported


    // choose saturated or unsaturated
    const double dM_b_dt = std::min(dM_b_dt_unsat, dM_b_dt_sat);

    const double dM_b = dM_b_dt * dt;


    const double x    = dM_b / c->cons_cold[DDD];
    const double dP   = x * c->cons_cold[SRR];
    // const double dE   = dM_b * c_s_cold * c_s_cold; // doesn't include kinetic energy
    const double dE   = x * c->cons_cold[TAU];        
    const double dM_Z = x * c->cons_cold[ZZZ];

    // ======== Verify before transferring mass / energy ========= //
    // #ifndef NDEBUG

    if (!std::isfinite(dM_b_dt_unsat))
    {
        printf("------ ERROR in subgrid_thermal_conduction()------- \n");
        printf("Non-finite unsaturated conduction rate! \n");
        printf("dM_b_dt_unsat = %e \n", dM_b_dt_unsat);
        printf("T_hot         = %e \n", T_hot);
        printf("T_hot**(5/2)  = %e \n", std::pow(T_hot, 5./2.));
        printf("dr            = %e \n", dr);
        assert( std::isfinite(dM_b_dt_unsat) );

    }

    if (!std::isfinite(dM_b_dt_sat))
    {
        printf("------ ERROR in subgrid_thermal_conduction()------- \n");
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
        printf("------ ERROR in subgrid_thermal_conduction()------- \n");
        printf("Negative mass transfered! \n");
        printf("dM_b = %e \n", dM_b);
        assert(dM_b > 0);
    }

    if( x < 0)
    {
        printf("------ ERROR in subgrid_thermal_conduction()------- \n");
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
        printf("------ ERROR in subgrid_thermal_conduction()------- \n");
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

        printf("------ ERROR in subgrid_thermal_conduction()------- \n");
        printf("No cold energy to begin with! \n");
        printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        assert(0 < c->cons_cold[TAU]);
    }   


    if (dE >= c->cons_cold[TAU])
    {
        const double T_cold = calc_T(c->prim_cold);

        printf("------ ERROR in subgrid_thermal_conduction()------- \n");
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

    c->cons_cold[SRR] -= dP;
    c->cons_hot[SRR]  += dP;

    c->cons_cold[TAU] -= dE;
    c->cons_hot[TAU]  += dE;

    c->cons_cold[ZZZ] -= dM_Z;
    c->cons_hot[ZZZ]  += dM_Z;


    // ======== Verify post-conditions ========= //
    if (c->multiphase)
    {
        verify_multiphase_conditions(c, "subgrid_thermal_conduction()", "post");

        // if (E_int_from_cons(c->cons_hot) < 0)
        // {
        //     printf("------ERROR in subgrid_thermal_conduction() postconditions------- \n");
        //     printf("hot gas *internal* energy less than 0\n");
        //     printf("E_int (hot)  = %e \n", E_int_from_cons(c->cons_hot));
        //     printf("c->cons_hot[TAU]   = %e \n", c->cons_hot[TAU]);
        //     printf("c->cons[TAU]       = %e \n", c->cons[TAU]);

        //     assert( E_int_from_cons(c->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(c->cons_cold) < 0)
        // {
        //     printf("------ERROR in subgrid_thermal_conduction() postconditions------- \n");
        //     printf("cold gas *internal* energy less than 0\n");
        //     printf("E_int (cold)  = %e \n", E_int_from_cons(c->cons_cold));
        //     printf("c->cons_cold[TAU]   = %e \n", c->cons_cold[TAU]);
        //     printf("c->cons[TAU]        = %e \n", c->cons[TAU]);

        //     assert( E_int_from_cons(c->cons_cold) > 0 );
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

void verify_multiphase_conditions( const struct cell * c ,
                                   const char * fn_name ,
                                   const char * pre_post_specifier ,
                                   const char * cell_specifier ,
                                   const double rel_tol )
{
    #ifndef NDEBUG

    if (c->cons_hot[DDD] < 0)
    {
        printf("------ERROR in %s %sconditions------- \n", fn_name, pre_post_specifier);
        printf("%shot gas mass less than 0\n", cell_specifier);
        printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
        printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
        printf("c->cons[DDD]      = %e \n", c->cons[DDD]);
        fflush(stdout);

        assert( c->cons_hot[DDD] > 0 );
    }


    if (c->cons_cold[DDD] < 0)
    {
        printf("------ERROR in %s %sconditions------- \n", fn_name, pre_post_specifier);
        printf("%scold gas mass less than 0\n", cell_specifier);
        printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
        printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
        printf("c->cons[DDD]      = %e \n", c->cons[DDD]);
        fflush(stdout);

        assert( c->cons_cold[DDD] > 0 );
    }


    if ( (c->cons_hot[DDD] / c->cons[DDD]) > (1+rel_tol))
    {
        printf("------ERROR in %s %sconditions------- \n", fn_name, pre_post_specifier);
        printf("%shot gas mass greater than total mass\n", cell_specifier);
        printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
        printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
        printf("c->cons[DDD]      = %e \n", c->cons[DDD]);
        fflush(stdout);

        assert( (c->cons_hot[DDD] / c->cons[DDD]) < (1+rel_tol) );
    }


    if ( (c->cons_cold[DDD] / c->cons[DDD]) > (1+rel_tol))
    {
        printf("------ERROR in %s %sconditions------- \n", fn_name, pre_post_specifier);
        printf("%scold gas mass greater than total mass\n", cell_specifier);
        printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
        printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
        printf("c->cons[DDD]      = %e \n", c->cons[DDD]);
        fflush(stdout);

        assert( (c->cons_cold[DDD] / c->cons[DDD]) < (1+rel_tol) );
    }


    if ( std::abs(1-( (c->cons_cold[DDD] + c->cons_hot[DDD])/c->cons[DDD])) > rel_tol)
    {
        printf("------ERROR in %s %sconditions------- \n", fn_name, pre_post_specifier);
        printf("%scold mass + hot mass =/= total mass\n", cell_specifier);
        printf("c->cons_cold[DDD] = %e \n", c->cons_cold[DDD]);
        printf("c->cons_hot[DDD]  = %e \n", c->cons_hot[DDD]);
        printf("c->cons[DDD]      = %e \n", c->cons[DDD]);
        fflush(stdout);

        assert(  std::abs(1-( (c->cons_cold[DDD] + c->cons_hot[DDD])/c->cons[DDD])) <= rel_tol);
    }


    if (c->cons_hot[TAU] < 0)
    {
        printf("------ERROR in %s %sconditions------- \n", fn_name, pre_post_specifier);
        printf("%shot gas energy less than 0\n", cell_specifier);
        printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
        printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        printf("c->cons[TAU]      = %e \n", c->cons[TAU]);
        fflush(stdout);

        assert( c->cons_hot[TAU] > 0 );
    }


    if (c->cons_cold[TAU] < 0)
    {
        printf("------ERROR in %s %sconditions------- \n", fn_name, pre_post_specifier);
        printf("%scold gas energy less than 0\n", cell_specifier);
        printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
        printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        printf("c->cons[TAU]      = %e \n", c->cons[TAU]);
        fflush(stdout);

        assert( c->cons_cold[TAU] > 0 );
    }


    if ( std::abs(1-( (c->cons_cold[TAU] + c->cons_hot[TAU])/c->cons[TAU])) > rel_tol)
    {
        printf("------ERROR in %s %sconditions------- \n", fn_name, pre_post_specifier);
        printf("cold energy + hot energy =/= total energy\n");
        printf("c->cons_cold[TAU] = %e \n", c->cons_cold[TAU]);
        printf("c->cons_hot[TAU]  = %e \n", c->cons_hot[TAU]);
        printf("c->cons[TAU]      = %e \n", c->cons[TAU]);
        fflush(stdout);

        assert(  std::abs(1-( (c->cons_cold[TAU] + c->cons_hot[TAU])/c->cons[TAU])) <= rel_tol);
    }

    if (! std::isfinite(c->cons_hot[DDD]))
    {
        assert(std::isfinite(c->cons_hot[DDD]));
    }

    #endif
}


void Picard_iteration( const double * values_old ,
                             double * values_guess ,
                             double * values_new , // will be overwritten
                       const int num_values , 
                       std::function<void(const double *, const double *, double *)> lambda_apply_linsolve,
                       const char * calling_fn_name ,
                       const double iteration_tolerance , // default arg
                       const int num_iterations_max , // default arg
                       const double alpha_relaxation  // default arg
                       )
    // ============================================= //
    //
    //  Wraps the iterative solving of an implicit equation using Picard
    //  iteration with a relaxation parameter
    //
    //  Inputs:
    //    - values_old, values_guess, values_new
    //        - arrays of the variable to be solved for (all doubles)
    //        - all have equal length, `num_cells`
    //          variables to be passed to `lambda_apply_linsolve`
    //        - `values_new` is purely an output, and can be garbage on input
    //    - num_values
    //        - length of `values_*`
    //    - lambda_apply_linsolve
    //        - function that takes `values_old`, `values_guess` and `values_new`
    //          and creates a new guess in the place of `values_new`
    //        - typically a closure around another function, made by a lambda                       const char * calling_fn_name , 
    //    - calling_fn_name 
    //        - a string to specify when you called Picard_iteration,
    //          for use within debugging/error print statements
    //    - iteration_tolerance
    //        - try to iterate until each new value changes by less than
    //          a factor of `iteration_tolerance` from the guess
    //          e.g. if `iteration_tolerance = 1.05`, you want a less than 5% 
    //          change between `guess` and `new`.
    //    - num_iterations_max
    //        - cutoff the iteration after `num_iterations_max` no matter what
    //        - ensures reasonable runtime, but indicates potentially large
    //          numberical errors
    //    - alpha_relaxation
    //        - the relaxation parameter, such that:
    //            - alpha = 1.0: just use new value returned by `lambda_apply_linsolve`
    //            - alpha = 0.0: never updates `values_guess`
    //                           (basically just re-runs )
    //      
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - affects `values_guess` and `values_new`
    //
    // ============================================= // 
{
    // ======== Verify pre-conditions ========= //
    if( num_values <= 0 )
    {
        printf("------ERROR in Picard_iteration() preconditions------- \n");
        printf("`num_values` must be greater than 0, not %d\n",
               num_values);
        assert( num_values > 0 );
    }

    if( iteration_tolerance <= 0 )
    {
        printf("------ERROR in Picard_iteration() preconditions------- \n");
        printf("`iteration_tolerance` must be greater than 0, not %e\n",
               iteration_tolerance);
        assert( iteration_tolerance > 0 );
    }

    if( num_iterations_max <= 0 )
    {
        printf("------ERROR in Picard_iteration() preconditions------- \n");
        printf("`num_iterations_max` must be greater than 0, not %d\n",
               num_iterations_max);
        assert( num_iterations_max > 0 );
    }

    if( (alpha_relaxation <= 0) || (alpha_relaxation >= 1))
    {
        printf("------ERROR in Picard_iteration() preconditions------- \n");
        printf("`alpha_relaxation` must be in (0, 1), not %f\n",
               alpha_relaxation);
        assert((alpha_relaxation > 0) && (alpha_relaxation < 1));
    }


    // ======== Primary Code ========= //

    for( int j = 0; j < num_iterations_max ; ++j)
    {
        // printf("\n=======implicit iteration j=%d ========\n\n", j);

        // given values_old and values_guess, overwrites values_new
        lambda_apply_linsolve(values_old, values_guess, values_new);

        double log_guess_old[num_values];
        double log_guess_new[num_values];

        bool stop_early = false;
        for( int i=0 ; i<num_values ; ++i )
        {
            log_guess_old[i] = std::log(values_guess[i]);
            log_guess_new[i] = std::log(values_new[i]);
        }   

        const double max_log_distance = max_norm(log_guess_old, 
                                                 log_guess_new,
                                                 num_values);

        if( max_log_distance < std::abs(std::log(iteration_tolerance)) )
        {
            stop_early = true;
            if(j > 5)
            {
                printf("num iterations = %d in `Picard_iteration` called from `%s`\n",
                    j, calling_fn_name);
                fflush(stdout);
            }
        }
                // printf("max_distance = %e in %s\n", std::exp(max_log_distance),
                    // calling_fn_name);
        if( j == (num_iterations_max-1) )
        {
            if( !stop_early )
            {

                throw ImplicitSolverFailedToConvergeError(calling_fn_name,
                                                          "Picard_iteration");
            }
        }

        for( int i=0 ; i<num_values ; ++i )
        {
            values_guess[i] =     alpha_relaxation  * values_new[i] 
                             + (1-alpha_relaxation) * values_guess[i];
        }

        if( stop_early )
        {
            break;
        }
    }

}


void Anderson_acceleration( const double * values_old ,
                                  double * values_guess ,
                                  double * values_new , // will be overwritten
                            const int num_values , 
                            std::function<void(const double *, const double *, double *)> lambda_apply_linsolve,
                            const char * calling_fn_name ,
                            const double iteration_tolerance , // default arg
                            const int num_iterations_max, // default arg
                            const int history_max
                           )
    // ============================================= //
    //
    //  Wraps the iterative solving of an implicit equation using Anderson-
    //  accelerated iteration.
    //
    //  Inputs:
    //    - values_old, values_guess, values_new
    //        - arrays of the variable to be solved for (all doubles)
    //        - all have equal length, `num_cells`
    //          variables to be passed to `lambda_apply_linsolve`
    //        - `values_new` is purely an output, and can be garbage on input
    //    - num_values
    //        - length of `values_*`
    //    - lambda_apply_linsolve
    //        - function that takes `values_old`, `values_guess` and `values_new`
    //          and creates a new guess in the place of `values_new`
    //        - typically a closure around another function, made by a lambda
    //    - calling_fn_name 
    //        - a string to specify when you called Anderson_acceleration,
    //          for use within debugging/error print statements
    //    - iteration_tolerance
    //        - try to iterate until each new value changes by less than
    //          a factor of `iteration_tolerance` from the guess
    //          e.g. if `iteration_tolerance = 1.05`, you want a less than 5% 
    //          change between `guess` and `new`.
    //    - num_iterations_max
    //        - cutoff the iteration after `num_iterations_max` no matter what
    //        - ensures reasonable runtime, but indicates potentially large
    //          numberical errors
    //    - history_max
    //        - use up to the most recent `history_max` for Anderson acceleration
    //
    //      
    //  Returns:
    //       void
    //
    //  Side effects:
    //    - affects `values_guess` and `values_new`
    //
    // ============================================= // 
{
    // ======== Verify pre-conditions ========= //
    if( num_values <= 0 )
    {
        printf("------ERROR in Anderson_acceleration() preconditions------- \n");
        printf("`num_values` must be greater than 0, not %d\n",
               num_values);
        assert( num_values > 0 );
    }

    if( iteration_tolerance <= 0 )
    {
        printf("------ERROR in Anderson_acceleration() preconditions------- \n");
        printf("`iteration_tolerance` must be greater than 0, not %e\n",
               iteration_tolerance);
        assert( iteration_tolerance > 0 );
    }

    if( num_iterations_max <= 0 )
    {
        printf("------ERROR in Anderson_acceleration() preconditions------- \n");
        printf("`num_iterations_max` must be greater than 0, not %d\n",
               num_iterations_max);
        assert( num_iterations_max > 0 );
    }

    if( history_max <= 0 )
    {
        printf("------ERROR in Anderson_acceleration() preconditions------- \n");
        printf("`history_max` must be greater than 0, not %d\n",
               history_max);
        assert( history_max > 0 );
    }

    // ======== Primary Code ========= //

    double *R = new double [history_max * num_values];
    double G[num_values][history_max];


    double xi[1][history_max];


    for( int j = 0; j < num_iterations_max ; ++j)
    {
        printf("\n=======anderson implicit iteration j=%d ========\n\n", j);
        int history;
        if( (j+1) < history_max )
        {
            // history is 1-index because it is a counting number
            history = j + 1; 
        }
        else
        {
            history = history_max;
        }

        // these must be reset *each time*
        double zeros[num_values][1];
        for( int i = 0 ; i<num_values ; ++i) zeros[i][0] = 0.0;

        double ones[1][history_max];
        for( int i = 0 ; i<history_max ; ++i) ones[0][i] = 1.0;

        double one[1][1] = {{1.0}};


        // given values_old and values_guess, overwrites values_new
        lambda_apply_linsolve(values_old, values_guess, values_new);

        // ============ quit early? ================= //
        // do this before the potentially-unneeded extra linear algebra

        double log_guess_old[num_values];
        double log_guess_new[num_values];

        for( int i=0 ; i<num_values ; ++i )
        {
            log_guess_old[i] = std::log(values_guess[i]);
            log_guess_new[i] = std::log(values_new[i]);
        }   

        const double max_log_distance = max_norm(log_guess_old, 
                                                 log_guess_new,
                                                 num_values);

        if( max_log_distance < std::abs(std::log(iteration_tolerance)) )
        {

            break; // stop early
        }
        else
        {
            if( j == (num_iterations_max-1) )
            {
                throw ImplicitSolverFailedToConvergeError(calling_fn_name,
                                                          "Anderson");
            }
        }

        // ============ linear algebra to generate a new guess ============= //


        // index that will cycle through columns of R, G using modulus wrapping
        const int current_column = j % history_max;

        for( int i=0 ; i<num_values ; ++i)
        {
            G[i][current_column] = values_new[i];

            // note, I'm using column major ordering, because that makes it 
            // simpler to pass to LAPACKE
            R[(num_values*current_column) + i] = (values_new[i] - values_guess[i]) / values_new[i];
        }

        // copy R to a temperature matrix, because A is likely to be overwritten
        // within the lapack call
        double *A = new double [history * num_values];
        for( int i=0; i<history*num_values ; ++i) A[i] = R[i];

        // pass to LAPACKe_dgglse
        // fortran documentation:  http://www.netlib.org/lapack/explore-html-3.6.1/d0/d85/dgglse_8f_a131c0fa85b2fb29d385ce87b199bf9aa.html
        //
        // Notation notes:
        // Mine   | theirs
        // ----------------
        //  A     |   A    
        //  ones  |   B
        //  zeros |   c
        //  one   |   d
        //  xi    |   x

        // but note that if j < history_max, we only want to use a subset
        // of R, G and not the entire thing!
        //
        // Fortunately, since I set up R using column major ordering,
        // I only have to change the value I set for N. Any later rows will
        // simply be ignored during the constrained linear solve.


        lapack_int M = num_values;
        lapack_int N = history;
        lapack_int P = 1;
        lapack_int LDA = M;
        lapack_int LDB = 1;

        lapack_int INFO; // result status: 0 for success

        INFO = LAPACKE_dgglse(
            LAPACK_COL_MAJOR,
            M, 
            N,
            P,
            A,
            LDA,
            *ones,
            LDB,
            *zeros,
            *one,
            *xi);

        if( INFO != 0 )
        {
            printf("INFO non-zero, indicating error. (INFO = %d)\n", INFO);
            assert(INFO==0);
        }

        double sum = 0;
        for( int i=0 ; i<history ; ++i) sum += xi[0][i];
        if( (sum < .95) || (sum > 1.05))
        {
            printf("sum not equal to 1. (sum = %f)\n", sum);
            for( int i=0 ; i<history ; ++i) printf("xi[%d] = %f\n", i, xi[0][i]);
            printf("INFO = %d\n", INFO);
            printf("`one` = %f\n", one[0][0]);

            assert((sum>.95) && (sum<1.05));
        }
        // printf("`one` = %f\n", one[0][0]);


        for( int i=0 ; i<num_values ; ++i )
        {
            // for picard iteration, this is where the "interesting" bit was
            // (the relaxation parameter)

            // now use the xi to get a weighted average of G rows
            // and then use this as our actual next guess
            values_new[i] = 0;
            for( int k=0 ; k<history ; ++k)
            {
                values_guess[i] += G[i][k] * xi[0][k];
            }
        }

        printf("weights: \n\n");
        for( int i=0 ; i<history ; ++i) printf("xi[%d] = %f\n", i, xi[0][i]);
        printf("\n");
        printf("sum = %f\n", sum);
        printf("\n\n\n");
        printf("current column = %d\n", current_column);

        delete[] A;

    }

    delete[] R;

}

double total_energy_of_all_cells(const struct domain * theDomain)
{
    const struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;

    double E_total = 0;
    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        const struct cell * c = &(theCells[i]);
        E_total += c->cons[TAU];
    }

    return E_total;
}

double total_energy_of_all_cells_from_prim(const struct domain * theDomain)
{
    const struct cell * theCells = theDomain->theCells;
    const int Nr = theDomain->Nr;
    const int Ng = theDomain->Ng;

    double E_total = 0;
    for( int i=Ng ; i<Nr-Ng ; ++i )
    {
        const struct cell * c = &(theCells[i]);
        const double mass = c->cons[DDD];
        const double v = c->prim[VRR];
        E_total += .5 * mass * v * v;
        E_total += mass * (1./(GAMMA_LAW-1.)) * c->prim[PPP] / c->prim[RHO]; 
    }

    return E_total;
}

EnergyChecker::EnergyChecker( const struct domain * theDomain,
                              const char * target_function_name )
    :   target_function_name(target_function_name)
{
    #ifndef NDEBUG
    this->set(theDomain);
    #endif
}

void EnergyChecker::set(const struct domain * theDomain)
{
    E_total_before = total_energy_of_all_cells(theDomain);
}

void EnergyChecker::check(const struct domain * theDomain)
{
    #ifndef NDEBUG
    const double E_total_after = total_energy_of_all_cells(theDomain);

    const double dE = E_total_after - E_total_before;

    const double rel_tol = 1e-10;
    const double abs_tol = 1e51 * rel_tol;

    if( std::abs(dE) > abs_tol )
    {
        printf("ERROR: dE exceeds abs_tol after '%s' \n",
            target_function_name);
        printf("E before: %e \n", E_total_before);
        printf("E after:  %e \n", E_total_after);
        printf("dE:       %e \n", dE);
        fflush(stdout);

        assert(std::abs(dE) < abs_tol);
    }

    if( std::abs((dE / E_total_before)) > rel_tol )
    {
        printf("ERROR: dE exceeds rel_tol after %s \n",
            target_function_name);
        printf("E before: %e \n", E_total_before);
        printf("E after:  %e \n", E_total_after);
        printf("dE:       %e \n", dE);
        printf("dE (%%):   %f \n", 1 - (dE / E_total_before));
        fflush(stdout);



        assert( std::abs( (dE / E_total_before)) < rel_tol );
    }
    #endif
}



