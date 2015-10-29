
#include <cmath>
#include <assert.h>
#include <stdlib.h>
extern "C" {
#include <grackle.h>
}

#include "../structure.H"
#include "../cooling.H"
#include "../geometry.H" // get_dV
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

    #ifndef NDEBUG
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
    #endif
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

double E_int_from_cons( const struct cell * c )
{

    const double mass  = c->cons[DDD];
    const double vr    = c->cons[SRR] / mass;
    const double E_int = c->cons[TAU] - .5*mass*vr*vr;

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
    #ifndef NDEBUG
    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        Ustar[q] = prim[q]*Ustar[DDD];
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

void source( const double * prim , double * cons , const double * grad , 
                   double * dE_cool , 
             const double rp , const double rm ,
             const double dV , const double dt ,
             const double R_shock, Cooling * cooling )
{
    const double P = prim[PPP];
    const double r = .5 * (rp + rm);

    // 1st order contribution
    cons[SRR] += 4 * M_PI * P * (std::pow(rp, 2.) - std::pow(rm, 2.)) * dt;
 
    // 2nd order contribution:
    cons[SRR] += 8*M_PI*grad[PPP]
                  *(     (std::pow(rp,3.) - std::pow(rm,3.))/3. 
                     + r*(std::pow(rp,2.) - std::pow(rm,2.))/2. ); 


    if( cooling->with_cooling == true )
    {
        bool cached_cooling = false;
        double dr = rp - rm;
        if( r > (R_shock+(10*dr)) ) cached_cooling = true;

        *dE_cool   = cooling->calc_cooling(prim, cons, dt , cached_cooling ) * dV;  
        cons[TAU] += *dE_cool; 
    }
    else
    {
        *dE_cool = 0;
    }
}

void conduction( struct cell * cL , struct cell * cR, 
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


