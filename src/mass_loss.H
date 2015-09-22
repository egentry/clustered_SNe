
#ifndef _MASS_LOSS_H_
#define _MASS_LOSS_H_


#include <string>
#include <vector>
#include "structure.H"


class Mass_Loss
{

public:
    Mass_Loss( const std::string name);

    virtual void add_mass_loss( struct domain * theDomain ,
                                const double dt ) const = 0;


    virtual double get_ejecta_mass(   const double M_initial ) const;
    virtual double get_ejecta_mass_Z( const double M_initial ) const;

    virtual double get_wind_mass( const double M_initial ) const;

    const std::string name;

private:

};



class No_Mass_Loss : public Mass_Loss
{

public:

    No_Mass_Loss();

    void add_mass_loss( struct domain * theDomain ,
                        const double dt ) const {};

    double get_ejecta_mass(   const double M_initial ) const ;
    double get_ejecta_mass_Z( const double M_initial ) const ;

    double get_wind_mass( const double M_initial ) const ;

    static const std::string class_name;


private:


};


class Uniform_Mass_Loss : public Mass_Loss
{

public:

    Uniform_Mass_Loss();

    void add_mass_loss( struct domain * theDomain , 
                        const double dt ) const;

    static const std::string class_name;

private:


};

class Disappear_Mass_Loss : public Mass_Loss
{

public:

    Disappear_Mass_Loss();

    void add_mass_loss( struct domain * theDomain ,
                        const double dt ) const {};


    double get_wind_mass( const double M_initial ) const ;

    static const std::string class_name;

private:


};

Mass_Loss * select_mass_loss( std::string mass_loss_name );

#endif

