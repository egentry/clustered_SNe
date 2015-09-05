
#ifndef _MASS_LOSS_H_
#define _MASS_LOSS_H_


#include <string>
#include <vector>
#include "structure.H"


class Mass_Loss
{

public:

    virtual void add_mass_loss( struct domain * theDomain ,
                                const double dt ) const = 0;


    virtual double get_ejecta_mass(   const double M_initial ) const;
    virtual double get_ejecta_mass_Z( const double M_initial ) const;

    virtual double get_wind_mass( const double M_initial ) const;

    virtual std::string get_name() const = 0;


};



class No_Mass_Loss : public Mass_Loss
{

public:

    void add_mass_loss( struct domain * theDomain ,
                        const double dt ) const {};

    double get_ejecta_mass(   const double M_initial ) const ;
    double get_ejecta_mass_Z( const double M_initial ) const ;

    double get_wind_mass( const double M_initial ) const ;

    std::string get_name() const ;



};


class Uniform_Mass_Loss : public Mass_Loss
{

public:

    void add_mass_loss( struct domain * theDomain , 
                        const double dt ) const;

    std::string get_name() const ;




};

class Disappear_Mass_Loss : public Mass_Loss
{

public:

    void add_mass_loss( struct domain * theDomain ,
                        const double dt ) const {};


    double get_wind_mass( const double M_initial ) const ;

    std::string get_name() const ;



};

Mass_Loss * select_mass_loss( std::string mass_loss_name );

#endif

