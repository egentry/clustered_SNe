
#include "structure.H"
#include "boundary.H"

void boundary( struct domain * theDomain ){

    struct cell * theCells = theDomain->theCells;
    int Nr = theDomain->Nr;


    struct cell * cB = &(theCells[0]);
    struct cell * cP = &(theCells[1]);
    cB->prim[RHO] = cP->prim[RHO];
    cB->prim[PPP] = cP->prim[PPP];
    cB->prim[ZZZ] = cP->prim[ZZZ];
    cB->prim[VRR] = 0;
    cB->wiph = 0;

    cB = &(theCells[Nr-1]);
    cP = &(theCells[Nr-2]);
    cB->prim[RHO] = cP->prim[RHO];
    cB->prim[PPP] = cP->prim[PPP];
    cB->prim[ZZZ] = cP->prim[ZZZ];
    cB->prim[VRR] = 0;
    cB->wiph = 0;

/*

    struct cell * cB = theCells[0];
    struct cell * cP = theCells[1];
    cB->prim[RHO] = cP->prim[RHO];
    cB->prim[PPP] = cP->prim[PPP];

*/
/*
    int ABSORB_R0 = theDomain->theParList.Absorb_BC;
    if( ABSORB_R0 )
    {
        struct cell * c3 = &(theCells[jk][1]);
        struct cell * c4 = &(theCells[jk][0]);
        for( int q=0 ; q<NUM_Q ; ++q )
        {
            c4->prim[q] = c3->prim[q];
        }    
    }    
*/

}
