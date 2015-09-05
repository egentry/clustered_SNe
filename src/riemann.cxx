
#include <cmath>
#include <assert.h>

#include "structure.H"
#include "misc.H" // prim2cons, cons2prim, 
#include "riemann.H"
#include "Hydro/euler.H" // flux, getUstar
#include "geometry.H" // get_dV

enum{_HLL_,_HLLC_};

static int riemann_solver = 0;

static double PRE_FLOOR = 0.0;


void setRiemannParams( const struct domain * theDomain )
{
    riemann_solver = theDomain->theParList.Riemann_Solver;
    PRE_FLOOR = theDomain->theParList.Pressure_Floor;
}


void riemann( struct cell * cL , struct cell * cR, 
              const double r , const double dA , const double dt )
{

    const double dAdt = dA * dt;
    double primL[NUM_Q];
    double primR[NUM_Q];

    const double drL = .5*cL->dr;
    const double drR = .5*cR->dr;

    for( int q=0 ; q<NUM_Q ; ++q )
    {
        primL[q] = cL->prim[q] + cL->grad[q]*drL;
        primR[q] = cR->prim[q] - cR->grad[q]*drR;
    }

    // ======== Verify pre-conditions ========= //
    #ifndef NDEBUG
    if (primL[PPP] < PRE_FLOOR)
    {
        printf("------ERROR in riemann()------- \n");
        printf("primL[%d] = %e \n", PPP, primL[PPP]);
        printf("Floor should be: %e \n", PRE_FLOOR);
        printf("cL->prim[%d] = %e \n", PPP, cL->prim[PPP]);
        printf("Occured at r = %e \n", r);
        printf("drL = %e \n", drL);
        printf("drR = %e \n", drR);
        printf("cL->grad[PPP] = %e \n", cL->grad[PPP]);

        assert(primL[PPP] > PRE_FLOOR);
    }

    if (primR[PPP] < PRE_FLOOR)
    {
        printf("------ERROR in riemann() ------- \n");
        printf("primR[%d] = %e \n", PPP, primR[PPP]);
        printf("Floor should be: %e \n", PRE_FLOOR);
        printf("cR->prim[%d] = %e \n", PPP, cR->prim[PPP]);
        printf("Occured at r = %e \n", r);
        printf("drL = %e \n", drL);
        printf("drR = %e \n", drR);
        printf("cR->grad[PPP] = %e \n", cR->grad[PPP]);

        assert(primR[PPP] > PRE_FLOOR);
    }
    #endif

    // continue with Riemann solver

    double Sl,Sr,Ss;

    vel( primL , primR , &Sl , &Sr , &Ss );

    double Fl[NUM_Q];
    double Fr[NUM_Q];
    double Ul[NUM_Q];
    double Ur[NUM_Q];

    double Flux[NUM_Q];

    const double w = cL->wiph;

    if( w < Sl ){
        flux( primL , Fl );
        prim2cons( primL , Ul , 1.0 );  // using dV = 1 gives the cons per unit volume

        for( int q=0 ; q<NUM_Q ; ++q ){
            Flux[q] = Fl[q] - w*Ul[q];

            #ifndef NDEBUG
            if( !std::isfinite(Flux[q]) )
            {
                printf("------ERROR in riemann()------- \n");
                printf("Bad flux in part 0 of riemann() \n");
                assert( std::isfinite(Flux[q]) );
            }
            #endif
        }
    }
    else if( w > Sr )
    {
        flux( primR , Fr );
        prim2cons( primR , Ur , 1.0 );

        for( int q=0 ; q<NUM_Q ; ++q )
        {
            Flux[q] = Fr[q] - w*Ur[q];
            #ifndef NDEBUG
            if( !std::isfinite(Flux[q]) )
            {
                printf("------ERROR in riemann()------- \n");
                printf("Bad flux in part 1 of riemann() \n");
                assert( std::isfinite(Flux[q]) );
            }
            #endif
        }
    }
    else
    {
        if( riemann_solver == _HLL_ )
        {
            const double aL =  Sr;
            const double aR = -Sl;

            prim2cons( primL , Ul , 1.0 );
            prim2cons( primR , Ur , 1.0 );
            flux( primL , Fl );
            flux( primR , Fr );

            for( int q=0 ; q<NUM_Q ; ++q )
            {
                const double Fstar = ( aL*Fl[q] + aR*Fr[q] + 
                                       aL*aR*(Ul[q] - Ur[q]) ) / (aL + aR);
                const double Ustar = ( aR*Ul[q] + aL*Ur[q] +
                                       Fl[q]    - Fr[q]      ) / (aL + aR);

                Flux[q] = Fstar - w*Ustar;
                #ifndef NDEBUG
                if( !std::isfinite(Flux[q]) )
                {
                    printf("------ERROR in riemann()------- \n");
                    printf("Bad flux in part 2 of riemann() \n");
                    assert( std::isfinite(Flux[q]) );
                }
                #endif
            }
        }
        else
        {
            double Ustar[NUM_Q];
            double Uk[NUM_Q];
            double Fk[NUM_Q];
            if( w < Ss )
            {
                prim2cons( primL , Uk , 1.0 );
                getUstar( primL , Ustar , Sl , Ss ); 
                flux( primL , Fk ); 

                for( int q=0 ; q<NUM_Q ; ++q )
                {
                    Flux[q] = Fk[q] + Sl*( Ustar[q] - Uk[q] ) - w*Ustar[q];
                    #ifndef NDEBUG
                    if( !std::isfinite(Flux[q]) )
                    {
                        printf("------ERROR in riemann()------- \n");
                        printf("Bad flux in part 3 of riemann() \n");
                        assert( std::isfinite(Flux[q]) );
                    }
                    #endif
                }    
            }
            else
            {
                prim2cons( primR , Uk , 1.0 );
                getUstar( primR , Ustar , Sr , Ss ); 
                flux( primR , Fk ); 

                for( int q=0 ; q<NUM_Q ; ++q )
                {
                    Flux[q] = Fk[q] + Sr*( Ustar[q] - Uk[q] ) - w*Ustar[q];
                    #ifndef NDEBUG
                    if( !std::isfinite(Flux[q]) )
                    {
                        printf("------ERROR in riemann()------- \n");
                        printf("Bad flux in part 4, Flux[%d]=%e \n", q, Flux[q]);
                        printf("Fk[%d] = %e \n", q, Fk[q]);
                        printf("Sr    = %e \n", Sr);
                        printf("Ustar[%d] = %e \n", q, Ustar[q]);
                        printf("w = %e \n", w);
                        assert( std::isfinite(Flux[q]) );
                    }
                    #endif
                } 
            } 
        }
    }

    #ifndef NDEBUG
    for ( int q=0 ; q<NUM_Q ; ++q)
    {
        if( !std::isfinite(Flux[q]) )
        {
            printf("------ERROR in riemann()------- \n");
            printf("Non-finite Flux[%d]= %e at r = %e \n", q, Flux[q], r);
            printf("cL->dr = %e \n", cL->dr);
            printf("cR->dr = %e \n", cR->dr);
            assert( std::isfinite(Flux[q]) );
        }
    }
    #endif


   for( int q=0 ; q<NUM_Q ; ++q )
   {
        cL->cons[q] -= Flux[q]*dAdt;
        cR->cons[q] += Flux[q]*dAdt;
   }


   // ======== Verify post-conditions ========= //
    #ifndef NDEBUG
        double primL_tmp[NUM_Q];
        double primR_tmp[NUM_Q];

        const double rp_L = cL->riph;
        const double rm_L = rp_L - cL->dr;
        const double dV_L = get_dV( rp_L , rm_L );
        cons2prim( cL->cons , primL_tmp , dV_L );  // using dV = 1 gives the cons per unit volume -- counteracts dV = 1 above
        const double rp_R = cR->riph;
        const double rm_R = rp_R - cR->dr;
        const double dV_R = get_dV( rp_R , rm_R );
        cons2prim( cR->cons , primR_tmp , dV_R);
        for( int q=0 ; q<NUM_Q ; ++q)
        {
            if( !std::isfinite(primL_tmp[q]) )
            {
                printf("------ERROR in riemann()------- \n");
                printf("primL[%d] not finite by end of riemann() \n", q);
                printf("primL[%d] = %e \n", q, primL_tmp[q]);
                assert( std::isfinite(primL_tmp[q]) );
            }
            if( !std::isfinite(primR_tmp[q]) )
            {
                printf("------ERROR in riemann()------- \n");
                printf("primR[%d] not finite by end of riemann() \n", q);
                printf("primR[%d] = %e \n", q, primR_tmp[q]);
                assert( std::isfinite(primR_tmp[q]) );
            }
        }

        if( primL_tmp[PPP] < PRE_FLOOR )
        {
            printf("------ERROR in riemann()------- \n");
            printf("while preparing to exit \n");
            printf("left prim[%d] = %e \n", PPP, primL_tmp[PPP]);
            printf("expected P > PRE_FLOOR \n");
            assert(primL_tmp[PPP] < PRE_FLOOR);
        }
        if( primR_tmp[PPP] < PRE_FLOOR )
        {
            printf("------ERROR in riemann()------- \n");
            printf("while preparing to exit \n");
            printf("right prim[%d] = %e \n", PPP, primR_tmp[PPP]);
            printf("expected P > PRE_FLOOR \n");
            assert(primR_tmp[PPP] < PRE_FLOOR);
        }
    #endif

}


