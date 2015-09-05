
#include <cmath>
#include "structure.H"

double minmod( const double a , const double b , const double c )
{
    double m = a;
    if( a*b < 0.0 )                 m = 0.0;
    if( std::abs(b) < std::abs(m) ) m = b;
    if( b*c < 0.0 )                 m = 0.0;
    if( std::abs(c) < std::abs(m) ) m = c;
    return m;
}

void plm( struct domain * theDomain )
{

    struct cell * theCells = theDomain->theCells;
    const int Nr  = theDomain->Nr;
    const int PLM = theDomain->theParList.PLM;
    for( int i=0 ; i<Nr ; ++i )
    {
        int im = i-1;
        int ip = i+1;
        if( i==0 ) im = 0;
        if( i==Nr-1 ) ip = Nr-1;
              struct cell * c  = &(theCells[i]);
        const struct cell * cL = &(theCells[im]);
        const struct cell * cR = &(theCells[ip]);
        const double drL = cL->dr;
        const double drC = c->dr;
        const double drR = cR->dr;
        for( int q=0 ; q<NUM_Q ; ++q )
        {
            const double pL = cL->prim[q];
            const double pC = c->prim[q];
            const double pR = cR->prim[q];
            double sL = pC - pL;
            sL /= .5*( drC + drL );
            double sR = pR - pC;
            sR /= .5*( drR + drC );
            double sC = pR - pL;
            sC /= .5*( drL + drR ) + drC;
            c->grad[q] = minmod( PLM*sL , sC , PLM*sR );
        }
    }
}

