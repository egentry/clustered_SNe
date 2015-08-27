
#include <math.h>
#include <assert.h>
#include <stdlib.h>
extern "C" {
#include <grackle.h>
}

#include "../structure.H"
#include "../cooling.H"
#include "euler.H"

static double GAMMA_LAW = 0.0;
static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;

void setHydroParams( struct domain * theDomain )
{
    GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
    RHO_FLOOR = theDomain->theParList.Density_Floor;
    PRE_FLOOR = theDomain->theParList.Pressure_Floor;
}

void prim2cons( double * prim , double * cons , double dV )
{
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vr  = prim[VRR];
    double v2 = vr*vr;
    double gam = GAMMA_LAW;
    double rhoe = Pp/(gam-1.);

    cons[DDD] = rho*dV;
    cons[SRR] = rho*vr*dV;
    cons[TAU] = (.5*rho*v2 + rhoe)*dV;

    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        cons[q] = cons[DDD]*prim[q];
    }
}

void cons2prim( double * cons , double * prim , double dV , bool verbose)
{

    // // E    :    total energy / unit volume
    // // e    : internal energy / unit mass
    // // rhoe : internal energy / unit volume
    // // Pp   : pressure (why the second 'p'?)
    double rho = cons[DDD]/dV;
    double Sr  = cons[SRR]/dV;
    double E   = cons[TAU]/dV;

    double vr = Sr/rho;
    double v2 = vr*vr;
    double rhoe = E - .5*rho*v2;
    double gam = GAMMA_LAW;
    double Pp = (gam-1.)*rhoe;

    #ifndef NDEBUG
    if( rho<RHO_FLOOR ) 
    {
        printf("------ ERROR in cons2prim() -- RHO_FLOOR ------- \n");
        printf("rho = %e \n", rho );
        printf("Expected rho > %e \n", RHO_FLOOR);
        printf("dV = %e \n", dV);
        rho=RHO_FLOOR; 

        assert(rho>RHO_FLOOR);
    }
    #endif
    if( Pp < PRE_FLOOR )
    {
        if (verbose)
        {
            printf("------ ERROR in cons2prim()------- \n");
            printf("pressure should be above pressure floor! \n");
            printf("pressure       = %e \n", Pp);
            printf("pressure floor = %e \n", PRE_FLOOR);
            printf("dV  = %e \n", dV);
            printf("rho = %e \n", rho);
            printf("vr  = %e \n", vr);
        }
        // assert(Pp > 0);

        Pp = PRE_FLOOR;
    } 

    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[VRR] = vr;

    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        prim[q] = cons[q]/cons[DDD];
    }

}

void getUstar( double * prim , double * Ustar , double Sk , double Ss )
{

    double rho  = prim[RHO];
    double vr   = prim[VRR];
    double Pp   = prim[PPP];
    double v2   = vr*vr;

    double gam  = GAMMA_LAW;

    double rhoe = Pp/(gam-1.);

    double rhostar =  rho*(Sk - vr)/(Sk - Ss);
    double Pstar   =   Pp*(Ss - vr)/(Sk - Ss);
    double Us      = rhoe*(Sk - vr)/(Sk - Ss);

    Ustar[DDD] = rhostar;
    Ustar[SRR] = rhostar*Ss;
    Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss*(Ss - vr) + Pstar;

    // ======== Verify post-conditions ========= //
    // require finite fluxes
    #ifndef NDEBUG
    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        Ustar[q] = prim[q]*Ustar[DDD];
        if( !isfinite(Ustar[q]) )
        {
            printf("Ustar[%d] = %e in 'getUstar()' \n", q, Ustar[q]);
            printf("prim[%d]  = %e in 'getUstar()' \n", q, prim[q]);
            printf("Sk        = %20.10le in 'getUstar()' \n", Sk);
            printf("Ss        = %20.10le in 'getUstar()' \n", Ss);
            assert( isfinite(Ustar[q]) );
        }
    }
    #endif

}

void flux( double * prim , double * flux )
{
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vr  = prim[VRR];
    double v2  = vr*vr;
    double gam = GAMMA_LAW;
    double rhoe = Pp/(gam-1.);

    flux[DDD] = rho*vr;
    flux[SRR] = rho*vr*vr + Pp;
    flux[TAU] = (.5*rho*v2 + rhoe + Pp)*vr;

    for( int q=ZZZ ; q<NUM_Q ; ++q )
    {
        flux[q] = flux[DDD]*prim[q];
    }
}

void source( double * prim , double * cons , double * grad , 
             double rp , double rm , double dV , double dt , 
             code_units cooling_units, int With_Cooling)
{
    double Pp  = prim[PPP];
    double r  = .5*(rp+rm);
    double r2 = (rp*rp+rm*rm+rp*rm)/3.;
    cons[SRR] += 2.*Pp*(r/r2)*dV*dt; // this is just the 1st order contribution
    // r / r2 is effectively (rp^2 - rm^2) / (rp^3 - rm^3) = d(surface area) / dV
    // so we're doing cons[SRR] += Pp * (4*pi*(rp^2 - rm^2)) * dt

    // includes 2nd order contribution:
    cons[SRR] += 8*M_PI*grad[PPP]
                  *(     (pow(rp,3.) - pow(rm,3.))/3. 
                     + r*(pow(rp,2.) - pow(rm,2.))/2. ); 


    if( With_Cooling == 1)
    {
        cons[TAU] += calc_cooling(prim, cons, dt, cooling_units) * dV;  
    }
}


void vel( double * prim1 , double * prim2 , 
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

    double gam  = GAMMA_LAW;

    double P1   = prim1[PPP];
    double rho1 = prim1[RHO];
    double vn1  = prim1[VRR];

    double cs1  = sqrt(fabs(gam*P1/rho1));

    double P2   = prim2[PPP];
    double rho2 = prim2[RHO];
    double vn2  = prim2[VRR];

    double cs2  = sqrt(fabs(gam*P2/rho2));

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
        assert(0);
    }
    #endif
  
}

double mindt( double * prim , double w , double r , double dr )
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

    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vr  = prim[VRR];
    double gam = GAMMA_LAW;

    double cs  = sqrt(fabs(gam*Pp/rho));

    double maxvr  = cs + fabs( vr - w );
    double dt     = dr/maxvr;

    return dt;
}


