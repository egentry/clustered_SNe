
#include "structure.H"
#include "boundary.H"

void boundary( struct domain * theDomain ){

    struct cell * theCells = theDomain->theCells;
    int Nr = theDomain->Nr;

    int rank = theDomain->rank;
    int size = theDomain->size;


    if( rank==0 )
    {    
        struct cell * cB = &(theCells[0]);
        struct cell * cP = &(theCells[1]);
        cB->prim[RHO] = cP->prim[RHO];
        cB->prim[PPP] = cP->prim[PPP];
        cB->prim[ZZZ] = cP->prim[ZZZ];
        cB->prim[VRR] = 0;
        cB->wiph = 0;
    }

    if( rank==size-1 )
    {
        struct cell * cB = &(theCells[Nr-1]);
        struct cell * cP = &(theCells[Nr-2]);
        cB->prim[RHO] = cP->prim[RHO];
        cB->prim[PPP] = cP->prim[PPP];
        cB->prim[ZZZ] = cP->prim[ZZZ];
        cB->prim[VRR] = 0;
        cB->wiph = 0;
    }

/*

    struct cell * cB = theCells[0];
    struct cell * cP = theCells[1];
    cB->prim[RHO] = cP->prim[RHO];
    cB->prim[PPP] = cP->prim[PPP];

*/

}
