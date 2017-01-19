
#include <cmath>
#include <assert.h>

#include "structure.H"
#include "misc.H" // prim2cons, cons2prim, 
#include "riemann.H"
#include "Hydro/euler.H" // flux, getUstar, E_int_from_cons
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
    const double rel_tol = 1e-4; // relative tolerance for float comparisons
    if (cL->multiphase)
    {
            // print out diagnostics before continuuing
            printf("\n - Riemann (left; begin) - \n");

            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);
            printf("\n");
            printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

        if (cL->cons_hot[DDD] < 0)
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[left] hot gas mass less than 0\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( cL->cons_hot[DDD] > 0 );
        }

        if (cL->cons_cold[DDD] < 0)
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[left] cold gas mass less than 0\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( cL->cons_cold[DDD] > 0 );
        }

        if ( (cL->cons_hot[DDD] / cL->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[left] hot gas mass greater than total mass\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( (cL->cons_hot[DDD] / cL->cons[DDD]) < (1+rel_tol) );
        }

        if ( (cL->cons_cold[DDD] / cL->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[left] cold gas mass greater than total mass\n");
            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

            assert( (cL->cons_cold[DDD] / cL->cons[DDD]) < (1+rel_tol) );
        }

        if ( std::abs(1-( (cL->cons_cold[DDD] + cL->cons_hot[DDD])/cL->cons[DDD])) > rel_tol)
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[left] cold mass + hot mass =/= total mass\n");
            printf("cL->cons_cold[DDD]                      = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons_hot[DDD]                       = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] + cL->cons_hot[DDD]  = %e \n", cL->cons_cold[DDD] + cL->cons_hot[DDD]);
            printf("cL->cons[DDD]                           = %e \n", cL->cons[DDD]);
            printf("relative error  = %e \n", 1 - ( (cL->cons_cold[DDD] + cL->cons_hot[DDD]) / cL->cons[DDD]));

            assert(  std::abs(1-( (cL->cons_cold[DDD] + cL->cons_hot[DDD])/cL->cons[DDD])) <= rel_tol);
        }



        if (cL->cons_hot[TAU] < 0)
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[left] hot gas energy less than 0\n");
            printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

            assert( cL->cons_hot[TAU] > 0 );
        }

        if (cL->cons_cold[TAU] < 0)
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[left] cold gas energy less than 0\n");
            printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

            assert( cL->cons_cold[TAU] > 0 );
        }

        // if (E_int_from_cons(cL->cons_hot) < 0)
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
        //     printf("[left] hot gas *internal* energy less than 0\n");
        //     printf("[left] E_int (hot)  = %e \n", E_int_from_cons(cL->cons_hot));
        //     printf("cL->cons_hot[TAU]   = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons[TAU]       = %e \n", cL->cons[TAU]);

        //     assert( E_int_from_cons(cL->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(cL->cons_cold) < 0)
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
        //     printf("[left] cold gas *internal* energy less than 0\n");
        //     printf("[left] E_int (cold)  = %e \n", E_int_from_cons(cL->cons_cold));
        //     printf("cL->cons_cold[TAU]   = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons[TAU]        = %e \n", cL->cons[TAU]);

        //     assert( E_int_from_cons(cL->cons_cold) > 0 );
        // }


        // if ( (cL->cons_hot[TAU] / cL->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
        //     printf("[left] hot gas energy greater than total energy\n");
        //     printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

        //     assert( (cL->cons_hot[TAU] / cL->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( (cL->cons_cold[TAU] / cL->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
        //     printf("[left] cold gas energy greater than total energy\n");
        //     printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
        //     printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
        //     printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

        //     assert( (cL->cons_cold[TAU] / cL->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( std::abs(1-( (cL->cons_cold[TAU] + cL->cons_hot[TAU])/cL->cons[TAU])) > rel_tol)
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
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
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[right] hot gas mass less than 0\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( cR->cons_hot[DDD] > 0 );
        }

        if (cR->cons_cold[DDD] < 0)
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[right] cold gas mass less than 0\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( cR->cons_cold[DDD] > 0 );
        }

        if ( (cR->cons_hot[DDD] / cR->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[right] hot gas mass greater than total mass\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( (cR->cons_hot[DDD] / cR->cons[DDD]) < (1+rel_tol) );
        }

        if ( (cR->cons_cold[DDD] / cR->cons[DDD]) > (1+rel_tol))
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[right] cold gas mass greater than total mass\n");
            printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);

            assert( (cR->cons_cold[DDD] / cR->cons[DDD]) < (1+rel_tol) );
        }

        if ( std::abs(1-( (cR->cons_cold[DDD] + cR->cons_hot[DDD])/cR->cons[DDD])) > rel_tol)
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[right] cold mass + hot mass =/= total mass\n");
            printf("cR->cons_cold[DDD]                      = %e \n", cR->cons_cold[DDD]);
            printf("cR->cons_hot[DDD]                       = %e \n", cR->cons_hot[DDD]);
            printf("cR->cons_cold[DDD] + cR->cons_hot[DDD]  = %e \n", cR->cons_cold[DDD] + cR->cons_hot[DDD]);
            printf("cR->cons[DDD]                           = %e \n", cR->cons[DDD]);
            printf("relative error  = %e \n", 1 - ( (cR->cons_cold[DDD] + cR->cons_hot[DDD]) / cR->cons[DDD]));

            assert(  std::abs(1-( (cR->cons_cold[DDD] + cR->cons_hot[DDD])/cR->cons[DDD])) <= rel_tol);
        }




        if (cR->cons_hot[TAU] < 0)
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[right] hot gas energy less than 0\n");
            printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
            printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
            printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

            assert( cR->cons_hot[TAU] > 0 );
        }

        if (cR->cons_cold[TAU] < 0)
        {
            printf("------ERROR in riemann() preconditions------- \n");
            printf("[right] cold gas energy less than 0\n");
            printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
            printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
            printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

            assert( cR->cons_cold[TAU] > 0 );
        }

        // if (E_int_from_cons(cR->cons_hot) < 0)
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
        //     printf("[right] hot gas *internal* energy less than 0\n");
        //     printf("[right] E_int (hot)  = %e \n", E_int_from_cons(cR->cons_hot));
        //     printf("cR->cons_hot[TAU]    = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons[TAU]        = %e \n", cR->cons[TAU]);

        //     assert( E_int_from_cons(cR->cons_hot) > 0 );
        // }

        // if (E_int_from_cons(cR->cons_cold) < 0)
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
        //     printf("[right] cold gas *internal* energy less than 0\n");
        //     printf("[right] E_int (cold)  = %e \n", E_int_from_cons(cR->cons_cold));
        //     printf("cR->cons_cold[TAU]    = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons[TAU]         = %e \n", cR->cons[TAU]);

        //     assert( E_int_from_cons(cR->cons_cold) > 0 );
        // }

        // if ( (cR->cons_hot[TAU] / cR->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
        //     printf("[right] hot gas mass greater than total mass\n");
        //     printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

        //     assert( (cR->cons_hot[TAU] / cR->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( (cR->cons_cold[TAU] / cR->cons[TAU]) > (1+rel_tol))
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
        //     printf("[right] cold gas mass greater than total mass\n");
        //     printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);

        //     assert( (cR->cons_cold[TAU] / cR->cons[TAU]) < (1+rel_tol) );
        // }

        // if ( std::abs(1-( (cR->cons_cold[TAU] + cR->cons_hot[TAU])/cR->cons[TAU])) > rel_tol)
        // {
        //     printf("------ERROR in riemann() preconditions------- \n");
        //     printf("[right] cold energy + hot energy =/= total energy\n");
        //     printf("cR->cons_cold[TAU]                      = %e \n", cR->cons_cold[TAU]);
        //     printf("cR->cons_hot[TAU]                       = %e \n", cR->cons_hot[TAU]);
        //     printf("cR->cons_cold[TAU] + cR->cons_hot[TAU]  = %e \n", cR->cons_cold[TAU] + cR->cons_hot[TAU]);
        //     printf("cR->cons[TAU]                           = %e \n", cR->cons[TAU]);
        //     printf("relative error  = %e \n", 1 - ( (cR->cons_cold[TAU] + cR->cons_hot[TAU]) / cR->cons[TAU]));

        //     assert(  std::abs(1-( (cR->cons_cold[TAU] + cR->cons_hot[TAU])/cR->cons[TAU])) <= rel_tol);
        // }
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
    //             printf("------ERROR in riemann()------- \n");
    //             printf("x_hot_L not within (0,1 + rel_tol) \n");
    //             printf("x_hot_L = %e \n", x_hot_L);
    //             printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
    //             printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
    //             printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

    //             assert( (x_hot_L > 0) && (x_hot_L < 1 + rel_tol) );
    //         }

    //         if( (y_hot_L > 1 + rel_tol)  || (y_hot_L < 0))
    //         {
    //             printf("------ERROR in riemann()------- \n");
    //             printf("y_hot_L not within (0,1 + rel_tol) \n");
    //             printf("y_hot_L = %e \n", y_hot_L);
    //             printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
    //             printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
    //             printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

    //             assert( (y_hot_L > 0) && (y_hot_L < 1 + rel_tol) );
    //         }

    //         if( (x_cold_L > 1 + rel_tol)  || (x_cold_L < 0))
    //         {
    //             printf("------ERROR in riemann()------- \n");
    //             printf("x_cold_L not within (0,1 + rel_tol) \n");
    //             printf("x_cold_L = %e \n", x_cold_L);
    //             printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
    //             printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
    //             printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

    //             assert( (x_cold_L > 0) && (x_cold_L < 1 + rel_tol) );
    //         }

    //         if( (y_cold_L > 1 + rel_tol)  || (y_cold_L < 0))
    //         {
    //             printf("------ERROR in riemann()------- \n");
    //             printf("y_cold_L not within (0,1 + rel_tol) \n");
    //             printf("y_cold_L = %e \n", y_cold_L);
    //             printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
    //             printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
    //             printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

    //             assert( (y_cold_L > 0) && (y_cold_L < 1 + rel_tol) );
    //         }


    //         if( std::abs(x_cold_L + x_hot_L - 1) >= rel_tol)
    //         {
    //             printf("------ERROR in riemann()------- \n");
    //             printf("x_cold_L + x_hot_L =/= 1 \n");
    //             printf("x_cold_L           = %e \n", x_cold_L);
    //             printf("x_hot_L            = %e \n", x_hot_L);
    //             printf("x_cold_L + x_hot_L = %e \n", x_cold_L + x_hot_L);

    //             assert( std::abs(x_cold_L + x_hot_L - 1) < rel_tol);
    //         }

    //         if( std::abs(y_cold_L + y_hot_L - 1) >= rel_tol)
    //         {
    //             printf("------ERROR in riemann()------- \n");
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
    //             printf("------ERROR in riemann()------- \n");
    //             printf("x_hot_R not within (0,1 + rel_tol) \n");
    //             printf("x_hot_R = %e \n", x_hot_R);
    //             printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
    //             printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
    //             printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

    //             assert( (x_hot_R > 0) && (x_hot_R < 1 + rel_tol) );
    //         }

    //         if( (y_hot_R > 1 + rel_tol)  || (y_hot_R < 0))
    //         {
    //             printf("------ERROR in riemann()------- \n");
    //             printf("y_hot_R not within (0,1 + rel_tol) \n");
    //             printf("y_hot_R = %e \n", y_hot_R);
    //             printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
    //             printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
    //             printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

    //             assert( (y_hot_R > 0) && (y_hot_R < 1 + rel_tol) );
    //         }

    //         if( (x_cold_R > 1 + rel_tol)  || (x_cold_R < 0))
    //         {
    //             printf("------ERROR in riemann()------- \n");
    //             printf("x_cold_R not within (0,1 + rel_tol) \n");
    //             printf("x_cold_R = %e \n", x_cold_R);
    //             printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
    //             printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
    //             printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);

    //             assert( (x_cold_R > 0) && (x_cold_R < 1 + rel_tol) );
    //         }

    //         if( (y_cold_R > 1 + rel_tol)  || (y_cold_R < 0))
    //         {
    //             printf("------ERROR in riemann()------- \n");
    //             printf("y_cold_R not within (0,1 + rel_tol) \n");
    //             printf("y_cold_R = %e \n", y_cold_R);
    //             printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
    //             printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
    //             printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

    //             assert( (y_cold_R > 0) && (y_cold_R < 1 + rel_tol) );
    //         }
            

    //         if( std::abs(x_cold_R + x_hot_R - 1) >= rel_tol)
    //         {
    //             printf("------ERROR in riemann()------- \n");
    //             printf("x_cold_R + x_hot_R =/= 1 \n");
    //             printf("x_cold_R           = %e \n", x_cold_R);
    //             printf("x_hot_R            = %e \n", x_hot_R);
    //             printf("x_cold_R + x_hot_R = %e \n", x_cold_R + x_hot_R);

    //             assert( std::abs(x_cold_R + x_hot_R - 1) < rel_tol);
    //         }

    //         if( std::abs(y_cold_R + y_hot_R - 1) >= rel_tol)
    //         {
    //             printf("------ERROR in riemann()------- \n");
    //             printf("y_cold_R + y_hot_R =/= 1 \n");
    //             printf("y_cold_R           = %e \n", y_cold_R);
    //             printf("y_hot_R            = %e \n", y_hot_R);
    //             printf("y_cold_R + y_hot_R = %e \n", y_cold_R + y_hot_R);

    //             assert( std::abs(y_cold_R + y_hot_R - 1) < rel_tol);
    //         }
    // }

    // const double E_int_L_before = E_int_from_cons(cL->cons);
    // const double E_int_R_before = E_int_from_cons(cR->cons);
    for( int q=0 ; q<NUM_Q ; ++q )
    {
        cL->cons[q] -= Flux[q] * dA * dt;
        cR->cons[q] += Flux[q] * dA * dt;
    }
    // const double E_int_L_after = E_int_from_cons(cL->cons);
    // const double E_int_R_after = E_int_from_cons(cR->cons);

    if (cL->multiphase)
    {

        printf("cL->x_hot     = %e \n", cL->x_hot);
        printf("cL->y_hot     = %e \n", cL->y_hot);
        printf("(1-cL->y_hot) = %e \n", (1-cL->y_hot));
        printf("cL->y_cold    = %e \n", cL->y_cold);

        printf("\n");

        // const double dM_L = -Flux[DDD] * dA * dt;
        // const double dP_L = -Flux[SRR] * dA * dt;
        // const double dE_L = -Flux[TAU] * dA * dt;

        // const double dE_int_L = E_int_L_after - E_int_L_before;
        // const double dE_kin_L = dE_L - dE_int_L;

        // const double z_hot_L  = cL->cons_hot[ZZZ]  / cL->cons_hot[DDD];
        // const double z_cold_L = cL->cons_cold[ZZZ] / cL->cons_cold[DDD];

        cL->cons_hot[DDD]  -= cL->x_hot  * Flux[DDD] * dA * dt;      
        cL->cons_cold[DDD] -= cL->x_cold * Flux[DDD] * dA * dt;      

        cL->cons_hot[SRR]  -= cL->x_hot  * Flux[SRR] * dA * dt;      
        cL->cons_cold[SRR] -= cL->x_cold * Flux[SRR] * dA * dt;    

        // cL->cons_hot[TAU]  += (y_hot_L  * dE_int_L) + (cL->x_hot  * dE_kin_L);      
        // cL->cons_cold[TAU] += (y_cold_L * dE_int_L) + (cL->x_cold * dE_kin_L); 

        // cL->cons_hot[TAU]  -= y_hot_L  * Flux[TAU] * dA * dt;      
        // cL->cons_cold[TAU] -= y_cold_L * Flux[TAU] * dA * dt; 

        cL->cons_hot[ZZZ]  -= cL->x_hot  * cL->z_hot  * Flux[DDD] * dA * dt; 
        cL->cons_cold[ZZZ] -= cL->x_cold * cL->z_cold * Flux[DDD] * dA * dt; 

        // printf("dE_int_L     = %e \n", dE_int_L);
        // printf("dE_kin_L     = %e \n", dE_kin_L);
        printf("cL->y_cold     = %e \n", cL->y_cold);
        printf("cL->y_hot      = %e \n", cL->y_hot);
        printf("cL->x_cold     = %e \n", cL->x_cold);
        printf("cL->x_hot      = %e \n", cL->x_hot);


        // printf("dE_int_L * cL->y_hot_L     = %e \n", dE_int_L * cL->y_hot);
        // printf("dE_int_L * cL->y_cold_L    = %e \n", dE_int_L * cL->y_cold);

        // printf("dE_kin_L * cL->x_hot_L     = %e \n", dE_kin_L * cL->x_hot);
        // printf("dE_kin_L * cL->x_cold_L    = %e \n", dE_kin_L * cL->x_cold);

        printf("\n");
    }


    if (cR->multiphase)
    {
        // const double dM_R = Flux[DDD] * dA * dt;
        // const double dP_R = Flux[SRR] * dA * dt;
        // const double dE_R = Flux[TAU] * dA * dt;

        // const double dE_int_R = E_int_R_after - E_int_R_before;
        // const double dE_kin_R = dE_R - dE_int_R;

        // const double z_hot_R  = cR->cons_hot[ZZZ]  / cR->cons_hot[DDD];
        // const double z_cold_R = cR->cons_cold[ZZZ] / cR->cons_cold[DDD];

        cR->cons_hot[DDD]  += cR->x_hot  * Flux[DDD] * dA * dt;      
        cR->cons_cold[DDD] += cR->x_cold * Flux[DDD] * dA * dt;      

        cR->cons_hot[SRR]  += cR->x_hot  * Flux[SRR] * dA * dt;      
        cR->cons_cold[SRR] += cR->x_cold * Flux[SRR] * dA * dt;    

        // cR->cons_hot[TAU]  += (y_hot_R  * dE_int_R) + (x_hot_R  * dE_kin_R);      
        // cR->cons_cold[TAU] += (y_cold_R * dE_int_R) + (x_cold_R * dE_kin_R);

        // cR->cons_hot[TAU]  += y_hot_R  * Flux[TAU] * dA * dt;      
        // cR->cons_cold[TAU] += y_cold_R * Flux[TAU] * dA * dt; 

        cR->cons_hot[ZZZ]  += cR->x_hot  * cR->z_hot  * Flux[DDD] * dA * dt;      
        cR->cons_cold[ZZZ] += cR->x_cold * cR->z_cold * Flux[DDD] * dA * dt; 
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

        if (cL->multiphase)
        {

            printf("\n - Riemann (left; end) - \n");

            printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
            printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
            printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);
            printf("\n");
            printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);

            if (cL->cons_hot[DDD] < 0)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[left] hot gas mass less than 0\n");
                printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
                printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
                printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert( cL->cons_hot[DDD] > 0 );
            }

            if (cL->cons_cold[DDD] < 0)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[left] cold gas mass less than 0\n");
                printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
                printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
                printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert( cL->cons_cold[DDD] > 0 );
            }

            if ( (cL->cons_hot[DDD] / cL->cons[DDD]) > (1+rel_tol))
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[left] hot gas mass greater than total mass\n");
                printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
                printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
                printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert( (cL->cons_hot[DDD] / cL->cons[DDD]) < (1+rel_tol) );
            }

            if ( (cL->cons_cold[DDD] / cL->cons[DDD]) > (1+rel_tol))
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[left] cold gas mass greater than total mass\n");
                printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
                printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
                printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert( (cL->cons_cold[DDD] / cL->cons[DDD]) < (1+rel_tol) );
            }

            if ( std::abs(1-( (cL->cons_cold[DDD] + cL->cons_hot[DDD])/cL->cons[DDD])) > rel_tol)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[left] cold mass + hot mass =/= total mass\n");
                printf("cL->cons_cold[DDD] = %e \n", cL->cons_cold[DDD]);
                printf("cL->cons_hot[DDD]  = %e \n", cL->cons_hot[DDD]);
                printf("cL->cons[DDD]      = %e \n", cL->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert(  std::abs(1-( (cL->cons_cold[DDD] + cL->cons_hot[DDD])/cL->cons[DDD])) <= rel_tol);
            }



            if (cL->cons_hot[TAU] < 0)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[left] hot gas energy less than 0\n");
                printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
                printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
                printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);
                printf("dE                 = %e \n", Flux[TAU] * dA * dt);

                assert( cL->cons_hot[TAU] > 0 );
            }

            if (cL->cons_cold[TAU] < 0)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[left] cold gas energy less than 0\n");
                printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
                printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
                printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);
                printf("dE                 = %e \n", Flux[TAU] * dA * dt);

                assert( cL->cons_cold[TAU] > 0 );
            }


            // if (E_int_from_cons(cL->cons_hot) < 0)
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[left] hot gas *internal* energy less than 0\n");
            //     printf("[left] E_int (hot)  = %e \n", E_int_from_cons(cL->cons_hot));
            //     printf("cL->cons_hot[TAU]   = %e \n", cL->cons_hot[TAU]);
            //     printf("cL->cons[TAU]       = %e \n", cL->cons[TAU]);

            //     assert( E_int_from_cons(cL->cons_hot) > 0 );
            // }

            // if (E_int_from_cons(cL->cons_cold) < 0)
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[left] cold gas *internal* energy less than 0\n");
            //     printf("[left] E_int (cold) = %e \n", E_int_from_cons(cL->cons_cold));
            //     printf("cL->cons_cold[TAU]  = %e \n", cL->cons_cold[TAU]);
            //     printf("cL->cons[TAU]       = %e \n", cL->cons[TAU]);

            //     assert( E_int_from_cons(cL->cons_cold) > 0 );
            // }



            // if ( (cL->cons_hot[TAU] / cL->cons[TAU]) > (1+rel_tol))
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[left] hot gas energy greater than total energy\n");
            //     printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            //     printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            //     printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);
            //     printf("dE                 = %e \n", Flux[TAU] * dA * dt);

            //     assert( (cL->cons_hot[TAU] / cL->cons[TAU]) < (1+rel_tol) );
            // }

            // if ( (cL->cons_cold[TAU] / cL->cons[TAU]) > (1+rel_tol))
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[left] cold gas energy greater than total energy\n");
            //     printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            //     printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            //     printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);
            //     printf("dE                 = %e \n", Flux[TAU] * dA * dt);

            //     assert( (cL->cons_cold[TAU] / cL->cons[TAU]) < (1+rel_tol) );
            // }


            // if ( std::abs(1-( (cL->cons_cold[TAU] + cL->cons_hot[TAU])/cL->cons[TAU])) > rel_tol)
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[left] cold energy + hot energy =/= total energy\n");
            //     printf("cL->cons_cold[TAU] = %e \n", cL->cons_cold[TAU]);
            //     printf("cL->cons_hot[TAU]  = %e \n", cL->cons_hot[TAU]);
            //     printf("cL->cons[TAU]      = %e \n", cL->cons[TAU]);
            //     printf("dE                 = %e \n", Flux[TAU] * dA * dt);

            //     assert(  std::abs(1-( (cL->cons_cold[TAU] + cL->cons_hot[TAU])/cL->cons[TAU])) <= rel_tol);
            // }

            printf(" - \n\n");

        }

        if (cR->multiphase)
        {
            if (cR->cons_hot[DDD] < 0)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[right] hot gas mass less than 0\n");
                printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
                printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
                printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert( cR->cons_hot[DDD] > 0 );
            }

            if (cR->cons_cold[DDD] < 0)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[right] cold gas mass less than 0\n");
                printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
                printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
                printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert( cR->cons_cold[DDD] > 0 );
            }

            if ( (cR->cons_hot[DDD] / cR->cons[DDD]) > (1+rel_tol))
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[right] hot gas mass greater than total mass\n");
                printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
                printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
                printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert( (cR->cons_hot[DDD] / cR->cons[DDD]) < (1+rel_tol) );
            }

            if ( (cR->cons_cold[DDD] / cR->cons[DDD]) > (1+rel_tol))
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[right] cold gas mass greater than total mass\n");
                printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
                printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
                printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert( (cR->cons_cold[DDD] / cR->cons[DDD]) < (1+rel_tol) );
            }

            if ( std::abs(1-( (cR->cons_cold[DDD] + cR->cons_hot[DDD])/cR->cons[DDD])) > rel_tol)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[right] cold mass + hot mass =/= total mass\n");
                printf("cR->cons_cold[DDD] = %e \n", cR->cons_cold[DDD]);
                printf("cR->cons_hot[DDD]  = %e \n", cR->cons_hot[DDD]);
                printf("cR->cons[DDD]      = %e \n", cR->cons[DDD]);
                printf("dM                 = %e \n", Flux[DDD] * dA * dt);

                assert(  std::abs(1-( (cR->cons_cold[DDD] + cR->cons_hot[DDD])/cR->cons[DDD])) <= rel_tol);
            }



            if (cR->cons_hot[TAU] < 0)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[right] hot gas energy less than 0\n");
                printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
                printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
                printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);
                printf("dE                 = %e \n", Flux[TAU] * dA * dt);

                assert( cR->cons_hot[TAU] > 0 );
            }

            if (cR->cons_cold[TAU] < 0)
            {
                printf("------ERROR in riemann() postconditions------- \n");
                printf("[right] cold gas energy less than 0\n");
                printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
                printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
                printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);
                printf("dE                 = %e \n", Flux[TAU] * dA * dt);

                assert( cR->cons_cold[TAU] > 0 );
            }

            // if (E_int_from_cons(cR->cons_hot) < 0)
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[right] hot gas *internal* energy less than 0\n");
            //     printf("[right] E_int (hot)  = %e \n", E_int_from_cons(cR->cons_hot));
            //     printf("cR->cons_hot[TAU]    = %e \n", cR->cons_hot[TAU]);
            //     printf("cR->cons[TAU]        = %e \n", cR->cons[TAU]);

            //     assert( E_int_from_cons(cR->cons_hot) > 0 );
            // }

            // if (E_int_from_cons(cR->cons_cold) < 0)
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[right] cold gas *internal* energy less than 0\n");
            //     printf("[right] E_int (cold)  = %e \n", E_int_from_cons(cR->cons_cold));
            //     printf("cR->cons_cold[TAU]    = %e \n", cR->cons_cold[TAU]);
            //     printf("cR->cons[TAU]         = %e \n", cR->cons[TAU]);

            //     assert( E_int_from_cons(cR->cons_cold) > 0 );
            // }

            // if ( (cR->cons_hot[TAU] / cR->cons[TAU]) > (1+rel_tol))
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[right] hot gas energy greater than total energy\n");
            //     printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
            //     printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
            //     printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);
            //     printf("dE                 = %e \n", Flux[TAU] * dA * dt);

            //     assert( (cR->cons_hot[TAU] / cR->cons[TAU]) < (1+rel_tol) );
            // }

            // if ( (cR->cons_cold[TAU] / cR->cons[TAU]) > (1+rel_tol))
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[right] cold gas energy greater than total energy\n");
            //     printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
            //     printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
            //     printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);
            //     printf("dE                 = %e \n", Flux[TAU] * dA * dt);

            //     assert( (cR->cons_cold[TAU] / cR->cons[TAU]) < (1+rel_tol) );
            // }


            // if ( std::abs(1-( (cR->cons_cold[TAU] + cR->cons_hot[TAU])/cR->cons[TAU])) > rel_tol)
            // {
            //     printf("------ERROR in riemann() postconditions------- \n");
            //     printf("[right] cold energy + hot energy =/= total energy\n");
            //     printf("cR->cons_cold[TAU] = %e \n", cR->cons_cold[TAU]);
            //     printf("cR->cons_hot[TAU]  = %e \n", cR->cons_hot[TAU]);
            //     printf("cR->cons[TAU]      = %e \n", cR->cons[TAU]);
            //     printf("dE                 = %e \n", Flux[TAU] * dA * dt);

            //     assert(  std::abs(1-( (cR->cons_cold[TAU] + cR->cons_hot[TAU])/cR->cons[TAU])) <= rel_tol);
            // }
        }

    #endif

}


