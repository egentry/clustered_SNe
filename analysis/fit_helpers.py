from __future__ import print_function, division

import numpy as np
from scipy import optimize

## Boilerplate path hack to give access to full clustered_SNe package
import sys, os
if __package__ is None:
    if os.pardir not in sys.path[0]:
        file_dir = os.getcwd()
        sys.path.insert(0, os.path.join(file_dir, 
                                        os.pardir, 
                                        os.pardir))

from clustered_SNe.analysis.constants import m_proton, pc, yr, M_solar, \
                                   metallicity_solar
from clustered_SNe.analysis.parse import Overview, RunSummary

from clustered_SNe.analysis.database_helpers import session, \
                                                Simulation, \
                                                Simulation_Status


#####################


def simplify_status(status):
    if status in {"Error", "Unphysical"}:
        return "bad"
    else:
        return "good"


class Aggregated_Results(object):
    """For use in fitting the momentum model"""
    def __init__(self):        

        metallicities   = []
        densities       = []
        masses          = []
        momenta         = []
        statuses        = []
        usable          = []
        num_SNe         = []

        for simulation in session.query(Simulation):
            status = session.query(Simulation_Status).get(simulation.id).status
            metallicities.append(simulation.metallicity)
            densities.append(simulation.background_density)
            masses.append(simulation.cluster_mass)
            momenta.append(simulation.momentum)
            statuses.append(simplify_status(status))
            num_SNe.append(simulation.num_SNe)
            
        self.metallicities   = np.array(metallicities)
        self.densities       = np.array(densities)
        self.masses          = np.array(masses)
        self.momenta         = np.array(momenta)
        self.statuses        = np.array(statuses)
        self.num_SNe         = np.array(num_SNe)

        self.usable = np.full_like(self.masses, False, dtype=bool)
        self.excluded = (  (self.densities > m_proton) \
                         & (self.masses < 10**3.1 * M_solar))

        self.momenta = self.momenta / self.masses / 1e5 # [km / s]
        self.masses = self.masses / M_solar # [M_solar]

        self.metallicities_1D = np.sort(list(set(self.metallicities)))
        self.densities_1D     = np.sort(list(set(self.densities)))
        self.masses_1D        = np.sort(list(set(self.masses)))

        # if one initialization is bad, all initializations are considered bad
        # (we can't average if at least one is unphysical)
        for metallicity in self.metallicities_1D:
            for density in self.densities_1D:
                for mass in self.masses_1D:
                    mask = np.isclose(self.metallicities, metallicity, atol=0) \
                         & np.isclose(self.densities,     density,     atol=0) \
                         & np.isclose(self.masses,        mass,        atol=0)
                    self.usable[mask] = np.all(self.statuses[mask] == "good")

        self.metallicities_3D, self.densities_3D, self.masses_3D = np.meshgrid(self.metallicities_1D,
                                                                               self.densities_1D,
                                                                               self.masses_1D, 
                                                                    indexing="ij") 

        self.momenta_3D = np.empty_like(self.masses_3D)

        for i, metallicity in enumerate(self.metallicities_1D):
            for j, density in enumerate(self.densities_1D):
                for k, mass in enumerate(self.masses_1D):

                    good_metallicities = (self.metallicities == metallicity)
                    good_densities     = (self.densities == density)
                    good_masses        = (self.masses == mass)
                    good =   good_metallicities \
                           & good_densities \
                           & good_masses \
                           & self.usable
                    
                    self.momenta_3D[i,j,k] = np.average(self.momenta[good])


        self.momenta_3D[self.momenta_3D == 0] = 1 # in case we take the log...

#####################


class Momentum_Model(object):
    
    def __init__(self, 
                 norm_low,            norm_high,
                 eta_metallicity_low, eta_metallicity_high, 
                 eta_density_low,     eta_density_high, 
                 eta_N_SNe_low,       eta_N_SNe_high):
        
        params_0 = (norm_low,            norm_high,
                    eta_metallicity_low, eta_metallicity_high,
                    eta_density_low,     eta_density_high,
                    eta_N_SNe_low,       eta_N_SNe_high)
        
        self.params_0 = params_0
        self.params   = params_0
    
    def __call__(self, metallicity, density, N_SNe):
        
        norm_low,            norm_high, \
        eta_metallicity_low, eta_metallicity_high, \
        eta_density_low,     eta_density_high, \
        eta_N_SNe_low,       eta_N_SNe_high = self.params
        
        x = (metallicity, density, N_SNe)
        
        return self.evaluate(x,
                             np.log10(norm_low),  np.log10(norm_high),
                             eta_metallicity_low, eta_metallicity_high,
                             eta_density_low,     eta_density_high,
                             eta_N_SNe_low,       eta_N_SNe_high)
        
    def __str__(self):
        s  = "p/(100 * M_sun * N_SNe) ~ min("
        s +=       "{0:.2e}"
        s +=       " * (Z / Z_sun)**{1:.2f}"
        s +=       " * (rho / m_p)**{2:.2f}"
        s +=       " * (N_SNe/1)**{3:.2f}"
        s += " , "
        s +=       "{4:.2e}"
        s +=       " * (Z / Z_sun)**{5:.2f}"
        s +=       " * (rho / m_p)**{6:.2f}"
        s +=       " * (N_SNe/1000)**{7:.2f}"
        s += ")"
        s = s.format(self.params[0], self.params[2], self.params[4], self.params[6],
                     self.params[1], self.params[3], self.params[5], self.params[7])
        
        return s
    def fit(self, x, y, fixed=None, **kwargs):
        """Fits the model using its initial self.params_0;
        Stores result as self.params
        

        Parameters
        --------
            x : locations of initial N data points; array: N x 3
                - metallicity 
                - density [g cm**-3]
                - number of SNe 
                - momentum efficiency [km / s]
            y : momentum efficiency [g cm**-3] at those N data points; array: Nx1
            fixed : mask of values to be held fixed; None or np.ndarray (dtype=bool)
                - fixed values set to True
                - free  values set to False
                - None for all values free
            **kwargs : to be passed to optimize.curve_fit; dict
                - does NOT get passed to the fitted function


        Returns
        -------
        popt : best-fit parameters; array: 8x1
        pcov : estimated parameter covariance; array: n_free x n_free
            - where: n_free = sum(~fixed)


        Side effects
        ------------
        Overwrites self.params with latest fit results
        
        Notes
        -----
        Fits self.evaluate
        
        covariance for normalizations is covariance in the *log*
        
        Doesn't return pcov entries corresponding to fixed values.
        This might cause errors if you forget to check if a parameter was fixed.
        In the future this should be fixed.
        """
        p0 = np.array(self.params_0)
        p0[0:2] = np.log10(p0[0:2])
        
        if fixed is None:
            fixed = np.full_like(p0, False, dtype=bool)
        elif isinstance(fixed, np.ndarray):
            if not fixed.dtype == bool:
                raise TypeError("fixed must be an np.ndarray of bools")
        else:
            raise TypeError("fixed must be type None or np.ndarray")
            
        def _f(x, *args):
            p0[~fixed] = args
            return self.evaluate(x, *p0)
        
        popt, pcov = optimize.curve_fit(_f, x, y, p0=p0[~fixed],
                                        xtol=.01,
                                        **kwargs)
        p0[~fixed] = popt
        popt = p0
        popt[0:2] = 10**popt[0:2] # we fit the log of the normalization
        self.params = tuple(popt)
        return popt, pcov

    
    @staticmethod
    def evaluate(x,
                log_norm_low,        log_norm_high,
                eta_metallicity_low, eta_metallicity_high,
                eta_density_low,     eta_density_high,
                eta_N_SNe_low,       eta_N_SNe_high):
        
        norm_low  = 10**log_norm_low
        norm_high = 10**log_norm_high
        
        metallicity, density, N_SNe = x
        
        # Normalize p_M_low at M=10**2 M_solar
        p_M_low  = norm_low  \
                    * (metallicity / metallicity_solar)**eta_metallicity_low \
                    * (density / m_proton)**eta_density_low \
                    * (N_SNe / (1))**eta_N_SNe_low
        
        # Normalize p_M_low at M=10**5 M_solar
        p_M_high = norm_high \
                    * (metallicity / metallicity_solar)**eta_metallicity_high \
                    * (density / m_proton)**eta_density_high \
                    * (N_SNe / (1000))**eta_N_SNe_high
        
        p_M = (p_M_low * p_M_high) / (p_M_low + p_M_high)
            
        return p_M

        




