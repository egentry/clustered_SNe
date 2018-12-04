
#ifndef _EULER_H_
#define _EULER_H_

#include <functional>

extern "C" {
#include <grackle.h>
}

#include "../cooling.H"

void setHydroParams( const struct domain * theDomain );
void prim2cons( const double * prim , double * cons , const double dV );
void cons2prim( const double * cons , double * prim , 
                const double dV , const bool verbose=true );
double E_kin_from_cons( const double * cons );
double E_int_from_cons( const double * cons );
double E_int_from_prim( const struct cell * c , 
                        const double gamma );
void getUstar( const double * prim , double * Ustar , 
               const double Sk , const double Ss );
void flux( const double * prim , double * flux );
void source( struct cell * c , 
             const double rp , const double rm , 
             const double dV , const double dt , 
             const double R_shock , Cooling * cooling );
void artificial_conduction( struct cell * cL , struct cell * cR, 
                            const double dA , const double dt );
void vel( const double * prim1 , const double * prim2 , 
          double * Sl , double * Sr , double * Ss );
double mindt( const double * prim , const double w , 
              const double r , const double dr );

double calc_T(const double * prim);

void calc_multiphase_prim(const struct cell * c, 
    double * prim_hot, double * prim_cold,
    double * V_hot,    double * V_cold );

void subgrid_thermal_conduction( struct cell * c , const double dt );

void verify_multiphase_conditions( const struct cell * c ,
                                   const char * fn_name ,
                                   const char * pre_post_specifier ,
                                   const char * cell_specifier = "" ,
                                   const double rel_tol = 1e-5 );

void thermal_conduction_implicit( struct domain * theDomain , const double dt );

void turbulent_diffusion_implicit( struct domain * theDomain ,
                                   const double dt );


void U_i_to_prim( const double * U , double * prim , const int i );

void thermal_conduction_apply_linsolve( const double * e_int_old ,
                                const double * e_int_guess ,
                                double * e_int_new ,
                                const double * r ,
                                const double * V ,
                                const double * rho ,
                                const double * mu ,
                                const double dt ,
                                const int num_cells ,
                                const int verbose );

[[deprecated]]
void spitzer_thermal_conduction_apply_linsolve( const double * e_int_old ,
                                const double * e_int_guess ,
                                double * e_int_new ,
                                const double * r ,
                                const double * V ,
                                const double * rho ,
                                const double * mu ,
                                const double dt ,
                                const int num_cells ,
                                const int verbose );

void turbulent_diffusion_apply_linsolve( const double * Y_old ,
                                         const double * Y_guess ,
                                         double * Y_new ,
                                         const double * r ,
                                         const double * Vol ,
                                         const double * v ,
                                         const double dt ,
                                         const double C_turbulent_diffusion ,
                                         const int num_cells ,
                                         const int verbose ,
                                         const bool allow_negatives = false );
double turbulent_diffusion_calculate_b(
    const double r_left ,
    const double r_right ,
    const double v_left ,
    const double v_right ,
    const double C_turbulent_diffusion 
     );

void turbulent_diffusion_get_fluxes( const double * densities , 
                                           double * mass_fluxes ,
                                     const double * r ,
                                     const double * v ,
                                     const double C_turbulent_diffusion ,
                                     const int num_cells ,
                                     const int verbose = 0 );


double max_norm( const double * A , const double * B, const int size );

void Picard_iteration( const double * values_old ,
                             double * values_guess ,
                             double * values_new ,
                       const int num_values ,
                       std::function<void(const double *, const double *, double*)> lambda_apply_linsolve,
                       const char * calling_fn_name , 
                       const double iteration_tolerance = 1. + 1e-3,
                       const int num_iterations_max = 100,
                       const double alpha_relaxation  = 0.5
                       );

void Anderson_acceleration( const double * values_old ,
                                  double * values_guess ,
                                  double * values_new , // will be overwritten
                            const int num_values , 
                            std::function<void(const double *, const double *, double *)> lambda_apply_linsolve,
                            const char * calling_fn_name ,
                            const double iteration_tolerance = 1. + 1e-3, // default arg
                            const int num_iterations_max = 100, // default arg
                            const int history_max = 8
                           );

double total_energy_of_all_cells(const struct domain * theDomain);
double total_energy_of_all_cells_from_prim(const struct domain * theDomain);

class EnergyChecker
{
    public:
        EnergyChecker( const struct domain * theDomain,
                       const char * target_function_name );
        void set( const struct domain * theDomain );
        void check( const struct domain * theDomain);

    protected:

    private:
        const char * target_function_name;
        double E_total_before;

};

class ImplicitSolverFailedToConvergeError : public std::exception
{
    private:
        std::string output_message;

    public:

    const char * implicit_solver_calling_fn_name;
    const char * iteration_type_name; // Picard vs. Anderson

    ImplicitSolverFailedToConvergeError(
      const char * implicit_solver_calling_fn_name,
      const char * iteration_type_name)
    :
      implicit_solver_calling_fn_name(implicit_solver_calling_fn_name),
      iteration_type_name(iteration_type_name)
    {
        char output[200];
        sprintf(output,
                "Output: Implicit solver (%s) failed to converge for `%s`",
                iteration_type_name,
                implicit_solver_calling_fn_name
                );
        output_message = std::string(output);
    }

    const char * what () const throw ()
    {

        return output_message.c_str();
    }
};

#endif
