
#include <iostream>
#include <string>

#include "../Hydro/euler.H" // prim2cons, cons2prim
#include "../boundary.H"
#include "../geometry.H" // get_moment_arm, get_dV
#include "../misc.H" // calc_dr

#include "initial_conditions.H"
#include "chevalier_ICs.H"
#include "cluster_SNe_ICs.H"
#include "ejecta_ICs.H"
#include "isentropic_ICs.H"
#include "messy_ICs.H"
#include "restart_ICs.H"
#include "shocktube_ICs.H"
#include "Thornton_parameter_study_ICs.H"
#include "uniform_ICs.H"

Initial_Conditions * select_initial_conditions( std::string IC_name )
{

    std::cout << "Using ICs: " << IC_name << std::endl;

    if ( IC_name.compare("chevalier") == 0 )
    {
        return new Chevalier_ICs;
    }
    if ( IC_name.compare("ejecta") == 0 )
    {
        return new Ejecta_ICs;
    }
    if ( IC_name.compare("cluster_SNe") == 0 )
    {
        return new Cluster_SNe_ICs;
    }
    if ( IC_name.compare("isentropic") == 0 )
    {
        return new Isentropic_ICs;
    }

    if ( IC_name.compare("messy") == 0 )
    {
        return new Messy_ICs;
    }

    if ( IC_name.compare("restart") == 0 )
    {
        return new Restart_ICs;
    }

    if ( IC_name.compare("shocktube") == 0 )
    {
        return new Shocktube_ICs;
    }

    if ( IC_name.compare("Thornton_parameter_study") == 0 )
    {
        return new Thornton_Parameter_Study_ICs;
    }

    if ( IC_name.compare("uniform") == 0 )
    {
        return new Uniform_ICs;
    }

    std::cerr << "Initial condition name didn't match known names" << std::endl;

}

void Initial_Conditions::setup_cells( struct domain * theDomain )
{

    int i;
    struct cell * theCells = theDomain->theCells;
    int Nr = theDomain->Nr;

    for( i=0 ; i<Nr ; ++i )
    {
        struct cell * c = &(theCells[i]);
        double rp = c->riph;
        double rm = rp - c->dr;
        c->wiph = 0.0; 
        double r = get_moment_arm( rp , rm );
        this->initial( c->prim , r ); 
        double dV = get_dV( rp , rm );
        prim2cons( c->prim , c->cons , dV );
        cons2prim( c->cons , c->prim , dV );
        c->P_old  = c->prim[PPP];
        c->dV_old = dV;
    }

    boundary( theDomain );

}

void Initial_Conditions::setup_grid( struct domain * theDomain )
{
    theDomain->Ng = NUM_G;
    int Num_R = theDomain->theParList.Num_R;
    int LogZoning = theDomain->theParList.LogZoning;

    double Rmin = theDomain->theParList.rmin;
    double Rmax = theDomain->theParList.rmax;

    int Nr = Num_R;

    theDomain->Nr = Nr;
    theDomain->theCells = (struct cell *) malloc( Nr*sizeof(struct cell));

    int i;

    double dx = 1./static_cast<double>(Num_R);
    double x0 = 0;
    double R0 = theDomain->theParList.LogRadius;
    for( i=0 ; i<Nr ; ++i )
    {
        double xp = x0 + (i+1.)*dx;
        double rp;
        if( LogZoning == 0 )
        {
            rp = Rmin + xp*(Rmax-Rmin);
        }else if( LogZoning == 1 )
        {
            rp = Rmin*pow(Rmax/Rmin,xp);
        }else
        {
            rp = R0*pow(Rmax/R0,xp) + Rmin-R0 + (R0-Rmin)*xp;
        }
        theDomain->theCells[i].riph = rp;
    }
    calc_dr( theDomain );

    this->setup_cells( theDomain );
};