from __future__ import print_function, division

import numpy as np
import pandas as pd
from scipy.integrate import quad

## Boilerplate path hack to give access to full SNe package
import sys, os
if __package__ is None:
    if os.pardir not in sys.path[0]:
        file_dir = os.path.dirname(__file__)
        sys.path.insert(0, os.path.join(file_dir, 
                                        os.pardir, 
                                        os.pardir, 
                                        os.pardir,))

from SNe.analysis.helper_functions import calculate_mean_molecular_weight, \
                                          calculate_entropy, \
                                          calculate_temperature

class SedovSolution(object):
    """docstring for SedovSolution"""
    def __init__(self, E_0, rho_0, metallicity, gamma=5/3):
        super(SedovSolution, self).__init__()

        self.E_0   = E_0   # blast energy
        self.rho_0 = rho_0 # background density (constant spatially) 

        self.mu = calculate_mean_molecular_weight(metallicity)

        self.gamma = gamma # adiabatic index
        self.j = 3         # spherical geometry
        self.w = 0         # uniform background

        if  self.w  >  ( self.j*(3-self.gamma) + 2*(self.gamma-1) ) / (self.gamma+1):
            raise ValueError("Vacuum case -- invalid choice of gamma")
        elif self.w == ( self.j*(3-self.gamma) + 2*(self.gamma-1) ) / (self.gamma+1):
            raise ValueError("Singular case -- invalid choice of gamma")

        # significant self-similatiry parameters, V
        self.V_0 = 2 / ( (self.j+2-self.w) * gamma )   # origin
        self.V_2 = 4 / ( (self.j+2-self.w)*(gamma+1) ) # shock front

        self.a = (self.j+2-self.w)*(gamma+1) / 4
        self.b = (gamma+1) / (gamma-1)
        self.c = (self.j+2-self.w) * gamma / 2
        self.d = (self.j+2-self.w) * (gamma+1) \
                / ( ((self.j+2-self.w)*(gamma+1)) - 2*(2+self.j*(gamma-1)))
        self.e = (2 + self.j*(gamma-1)) / 2

        self.alpha_0 = 2 / (self.j+2-self.w)
        self.alpha_2 = - (gamma-1) / ( 2*(gamma-1) + self.j - gamma*self.w )
        self.alpha_1 = ( (self.j+2-self.w)*gamma / (2 + self.j*(gamma-1)) ) \
                  * ( ( 2*(self.j*(2-gamma) - self.w)\
                     / (gamma*(self.j+2-self.w))**2 )     - self.alpha_2 )
        self.alpha_3 = (self.j-self.w) / ( 2*(gamma-1) + self.j - gamma*self.w )
        self.alpha_4 = self.alpha_1 * (self.j+2-self.w)*(self.j-self.w) \
                        / ( self.j*(2-gamma) - self.w )
        self.alpha_5 = ( self.w*(1+gamma) - 2*self.j ) \
                        / ( self.j*(2-gamma) - self.w )

        self.xi_0    = self.get_xi_0()

        # amortized results:
        self.E_kin = self.get_E_kin()
        self.E_int = self.get_E_int()
        self.E_tot = self.E_kin + self.E_int

        self.momentum_dimensionless = self.get_momentum_dimensionless()

    def x_1(self,V):
        return self.a * V
    def x_2(self,V):
        return self.b * (self.c*V - 1)
    def x_3(self,V):
        return self.d * (1 - self.e*V)
    def x_4(self,V):
        return self.b * (1 - (self.c/self.gamma)*V)

    def l(self, V):
        """Calculates l := r / r_2 (where the subscript 2 denotes post-shock)
        """
        l = self.x_1(V)**(-self.alpha_0) \
          * self.x_2(V)**(-self.alpha_2) \
          * self.x_3(V)**(-self.alpha_1)
        return l

    def f(self, V):
        """Calculates f := v / v_2 (where the subscript 2 denotes post-shock)
        """
        f = self.x_1(V) * self.l(V)
        return f

    def g(self, V):
        """Calculates g := rho / rho_2 (where the subscript 2 denotes post-shock)
        """
        g = self.x_1(V)**(self.alpha_0*self.w) \
          * self.x_2(V)**(self.alpha_3 + self.alpha_2*self.w) \
          * self.x_3(V)**(self.alpha_4 + self.alpha_1*self.w) \
          * self.x_4(V)**(self.alpha_5)
        return g

    def h(self, V):
        """Calculates h := P / P_2 (where the subscript 2 denotes post-shock)
        """
        h = self.x_1(V)**(self.alpha_0*self.j) \
          * self.x_3(V)**(self.alpha_4 + self.alpha_1*(self.w-2)) \
          * self.x_4(V)**(1 + self.alpha_5)
        return h

    def dr_dV(self, V):
        """Calculate the jacobian factor dr/dV 
        """
        dr_dV = self.l(V) \
          * ( ( self.a       )*(-self.alpha_0)/self.x_1(V) \
            + ( self.b*self.c)*(-self.alpha_2)/self.x_2(V) \
            + (-self.d*self.e)*(-self.alpha_1)/self.x_3(V) )
        return dr_dV

    def get_momentum_dimensionless(self, verbose=False):
        """Calculate the dimensionless momentum
            
        Follows the approach of http://cococubed.asu.edu/papers/kamm_2000.pdf

        Parameters
        ----------
            gamma : Optional[float]
                Adiabatic index 
            verbose: Optional[float]
                Print the integral results?

        Returns
        -------
            momentum : float
                dimensionless momentum 

        Side effects
        ------------
            None

        Notes
        -----
            All values are to be re-dimensionalized using their post-shock values
            i.e. for post-shock values denoted by the subscript '2':
                Momentum *= r_2**3 * rho_2 * v_2

        """

        dMomentum = lambda V: (4 * np.pi * self.l(V)**2 * self.dr_dV(V)) \
                    * self.g(V) * self.f(V)

        momentum  = quad(dMomentum, self.V_0, self.V_2)

        if verbose is True:
            print("momentum = {0:e} +/- {1:e}".format(momentum[0], momentum[1]))

        return momentum[0]

    def get_momentum(self, time=None):
        if time is None:
            try:
                time = self.time
            except AttributeError:
                print("self.time must be set [e.g. through update()] or passed as an option")
                raise

        r_shock = self.get_r_shock(time=time)
        rho_2   = self.get_rho_2()
        u_2     = self.get_u_2(time=time)

        momentum = self.momentum_dimensionless * r_shock**3 * rho_2 * u_2

        return momentum


    def get_E_kin_dimensionless(self, verbose=False):
        """Calculate the dimensionless kinetic energy
            
        Follows the approach of http://cococubed.asu.edu/papers/kamm_2000.pdf

        Parameters
        ----------
            gamma : Optional[float]
                Adiabatic index 
            verbose: Optional[float]
                Print the integral results?


        Returns
        -------
            E_kin : float
                dimensionless kinetic energy


        Side effects
        ------------
            None


        Notes
        -----
            All values are to be re-dimensionalized using their post-shock values
            i.e. for post-shock values denoted by the subscript '2':
                E_kin    *= r_2**3 * rho_2 * v_2**2

        """

        dE_kin = lambda V: (4 * np.pi * self.l(V)**2 * self.dr_dV(V)) \
                * 0.5 * self.g(V) * self.f(V)**2

        E_kin  = quad(dE_kin, self.V_0, self.V_2)

        if verbose is True:
            print("E_kin    = {0:e} +/- {1:e}".format(E_kin[0], E_kin[1]))

        return E_kin[0]


    def get_E_kin(self):
        time = 1
        r_shock = self.get_r_shock(time=time)
        rho_2   = self.get_rho_2()
        u_2     = self.get_u_2(time=time)

        E_kin = self.get_E_kin_dimensionless() * r_shock**3 * rho_2 * u_2**2
        return E_kin


    def get_E_int_dimensionless(self, verbose=False):
        """Calculate the dimensionless internal energy
            
        Follows the approach of http://cococubed.asu.edu/papers/kamm_2000.pdf

        Parameters
        ----------
            gamma : Optional[float]
                Adiabatic index 
            verbose: Optional[float]
                Print the integral results?

        Returns
        -------
            E_int : (float, float)
                dimensionless internal energy


        Side effects
        ------------
            None


        Notes
        -----
            All values are to be re-dimensionalized using their post-shock values
            i.e. for post-shock values denoted by the subscript '2':
                E_int    *= r_2**3 * P_2

        """

        dE_int = lambda V: (4 * np.pi * self.l(V)**2 * self.dr_dV(V)) \
                * self.h(V) / (self.gamma-1)

        E_int  = quad(dE_int, self.V_0, self.V_2)

        if verbose is True:
            print("E_int    = {0:e} +/- {1:e}".format(E_int[0], E_int[1]))

        return E_int[0]


    def get_E_int(self):
        time = 1
        r_shock = self.get_r_shock(time=time)
        P_2     = self.get_P_2(time=time)

        E_int = self.get_E_int_dimensionless() * r_shock**3 * P_2
        return E_int


    def get_r_shock(self, xi_0=None, time=None):
        """Get the radius of the shock front"""
        if xi_0 is None:
            try:
                xi_0 = self.xi_0
            except AttributeError:
                print("self.xi_0 must already be set, or passed as an option")
                raise

        if time is None:
            try:
                time = self.time
            except AttributeError:
                print("self.time must be set [e.g. through update()] or passed as an option")
                raise

        r_shock = xi_0 * (self.E_0 * time**2 / self.rho_0)**.2
        return r_shock  


    def get_u_shock(self, xi_0=None, time=None):
        """Get the velocity of the shock front (not fluid velocity!)"""
        if xi_0 is None:
            try:
                xi_0 = self.xi_0
            except AttributeError:
                print("self.xi_0 must already be set, or passed as an option")
                raise

        if time is None:
            try:
                time = self.time
            except AttributeError:
                print("self.time must be set [e.g. through update()] or passed as an option")
                raise

        u_shock = .4 * xi_0 * (self.E_0 / (self.rho_0 * time**3) )**.2
        return u_shock    
        
    def get_rho_2(self):
        """Get the post-shock density"""
        rho_2 = (self.gamma+1)/(self.gamma-1) * self.rho_0
        return rho_2

    def get_u_2(self, xi_0=None, time=None):
        """Get post-shock fluid velocity"""
        if xi_0 is None:
            try:
                xi_0 = self.xi_0
            except AttributeError:
                print("self.xi_0 must already be set, or passed as an option")
                raise

        if time is None:
            try:
                time = self.time
            except AttributeError:
                print("self.time must be set [e.g. through update()] or passed as an option")
                raise

        u_shock = self.get_u_shock(xi_0=xi_0, time=time)
        u_2 = (2 / (self.gamma+1)) * u_shock
        return u_2

    def get_P_2(self, xi_0=None, time=None):
        """Get post-shock fluid pressure"""
        if xi_0 is None:
            try:
                xi_0 = self.xi_0
            except AttributeError:
                print("self.xi_0 must already be set, or passed as an option")
                raise

        if time is None:
            try:
                time = self.time
            except AttributeError:
                print("self.time must be set [e.g. through update()] or passed as an option")
                raise

        u_shock = self.get_u_shock(xi_0=xi_0, time=time)

        P_2 = (2 / (self.gamma+1)) * self.rho_0 * u_shock**2
        return P_2

    def get_xi_0(self):
        """Uses the energy integral to the get the numerical factor xi_0

        xi_0 is the factor defined in the terminology of 
        http://www.astronomy.ohio-state.edu/~ryden/ast825/ch5-6.pdf

        xi_0 ~ 1.17, but that's not good enough for precision work

        Side effects
        ------------
            None
        """

        E_kin_dimensionless = self.get_E_kin_dimensionless()
        E_int_dimensionless = self.get_E_int_dimensionless()

        # set the current time to an arbitrary time
        time    = 1e7

        # dimensionalize as if xi_0 = 0
        xi_0 = 1


        r_shock = self.get_r_shock(xi_0=xi_0, time=time)
        rho_2 = self.get_rho_2()
        u_2   = self.get_u_2(xi_0=xi_0, time=time)
        P_2   = self.get_P_2(xi_0=xi_0, time=time)

        E_kin = E_kin_dimensionless * r_shock**3 * rho_2 * u_2**2
        E_int = E_int_dimensionless * r_shock**3 * P_2

        xi_0 = (self.E_0 / (E_kin + E_int))**.2
        return xi_0

    def update(self, time):
        try:
            time_old = self.time
            reuse_solution = True
        except:
            # if there isn't an old time, we need to make a new solution
            reuse_solution = False

        self.time = time

        self.momentum = self.get_momentum()

        if reuse_solution:
            self.solution["radius"]      *= (time / time_old)**+0.4
            self.solution["velocity"]    *= (time / time_old)**-0.6
            self.solution["density"]     *= (time / time_old)**0
            self.solution["pressure"]    *= (time / time_old)**-1.2

            self.solution["temperature"] *= (time / time_old)**-1.2
            self.solution["c_ad"]        *= (time / time_old)**-0.6
            self.solution["energy"]      *= (time / time_old)**-1.2

            self.solution["m_int"]       *= (time / time_old)**+1.2
        else:
            solution = self.generate_solution(num=50000)

            self.solution                = solution
            self.solution["temperature"] = self.get_temperature()
            self.solution["c_ad"]        = np.sqrt(self.gamma \
                                            * self.solution["pressure"] \
                                            / self.solution["density"])

            self.solution["energy"]      = (  self.solution["pressure"] \
                                            / self.solution["density"] ) \
                                            / (self.gamma-1)

            radius = self.solution["radius"]
            radius = np.append(0, radius.values)
            dV = np.diff( (radius**3).cumsum() )
            self.solution["m_int"] = dV * self.solution["density"]


        self.solution["entropy"]     = self.get_entropy()

        return

    def get_temperature(self):
        temperature = calculate_temperature( self.solution["pressure"], 
                                             self.solution["density"], 
                                             self.mu)
        return temperature

    def get_entropy(self):
        """Calculates the current entropy

        Warning
        -------
            Expects that the temperature is already updated,
            but that is not explicitly enforced
        """
        entropy = calculate_entropy(self.solution["temperature"],
                                    self.solution["density"],
                                    self.mu)
        return entropy


    def generate_solution(self, num=100):
        """Generate the Sedov solution

        Builds the non-dimensional solution using:
            http://cococubed.asu.edu/papers/kamm_2000.pdf

        then dimensionalizes it using:
            http://www.astronomy.ohio-state.edu/~ryden/ast825/ch5-6.pdf

        Parameters
        ----------
            time : float
                time since explosion
            num : Optional[int]
                number of gridpoints to generate

        Returns
        -------
            radius : np.array[float]
                radius of grid points
                Not uniform!
            velocity : np.array[float]
                velocities of fluid (not shock velocity)
            density : np.array[float]
            pressure : np.array[float]
            temperature : np.array[float]
            entropy : np.array[float]


        Notes
        -----
        We build the solution using linear spacing on the similarity variable V,
        not on a uniform grid of radius. This results in a poor sampling at 
        the inner radii.


        If you need a uniform radial grid either interpolate,
        or better yet, generate a solution using sedov3.f from cococubed
        """
        # V           = np.linspace(self.V_0, self.V_2, num=num)
        V           = self.V_0 \
                        + (self.V_2 - self.V_0)*np.linspace(0,1,num=num)**7

        V           = V[1:]    # ignore the innermost location, where r=0
        radius      = self.l(V)
        velocity    = self.f(V)
        density     = self.g(V)
        pressure    = self.h(V)


        r_shock = self.get_r_shock()
        rho_2   = self.get_rho_2()
        u_2     = self.get_u_2()
        P_2     = self.get_P_2()

        radius   *= r_shock
        velocity *= u_2
        density  *= rho_2
        pressure *= P_2

        solution_array = np.array([radius, velocity, density, pressure]).transpose()
        df = pd.DataFrame.from_records(solution_array,
            columns=["radius", "velocity", "density", "pressure"])

        return df




# sedov = SedovSolution(1e51, 1e-24, .02)
# print("E_kin: ", sedov.E_kin)
# print("E_int: ", sedov.E_int)

# sedov.update(1e10)



