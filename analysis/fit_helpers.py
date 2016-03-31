from __future__ import print_function, division

from matplotlib import pyplot as plt


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
from clustered_SNe.analysis.parse import Overview, RunSummary, \
                                         parse_into_scientific_notation

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
        ids             = []

        for simulation in session.query(Simulation):
            status = session.query(Simulation_Status).get(simulation.id).status
            metallicities.append(simulation.metallicity)
            densities.append(simulation.background_density)
            masses.append(simulation.cluster_mass)
            momenta.append(simulation.momentum)
            statuses.append(simplify_status(status))
            num_SNe.append(simulation.num_SNe)
            ids.append(simulation.id)
            
        self.metallicities   = np.array(metallicities)
        self.densities       = np.array(densities)
        self.masses          = np.array(masses)
        self.momenta         = np.array(momenta)
        self.statuses        = np.array(statuses)
        self.num_SNe         = np.array(num_SNe)
        self.ids             = np.array(ids)

        self.usable = np.full_like(self.masses, False, dtype=bool)
        self.excluded = (  (self.densities > m_proton) \
                         & (self.masses < 10**3.1 * M_solar))

        self.momenta = self.momenta / 1e5 # [g km / s]
        self.masses  = self.masses  / M_solar # [M_solar]

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
                    
                    self.momenta_3D[i,j,k] = np.average(self.momenta[good] / self.masses[good])


        self.momenta_3D[self.momenta_3D == 0] = 1 # in case we take the log...

        self.add_overall_fit()

    def add_overall_fit(self, verbose=False):
        if verbose:
            print_device = sys.stdout
        else:
            print_device = open(os.devnull, mode="w")

        self.model_overall = Momentum_Model(1e3, 1e4,
                                            0, 0, 
                                           .25, .2,
                                           2, -.08)

        fixed = np.array([False, False, 
                          False, False, 
                          False, False, 
                          False, False])

        mask = self.usable & (self.momenta>0)
        xs = (self.metallicities[mask],
              self.densities[mask],
              self.num_SNe[mask])
        ys = self.momenta[mask] / (self.num_SNe[mask] * 100 * M_solar)


        popt, pcov = self.model_overall.fit(xs, ys, fixed=fixed)
        print("Overall Model:", file=print_device)
        print("params_0: ", self.model_overall.params_0, file=print_device)
        print("params:   ", self.model_overall.params, file=print_device)
        print(file=print_device)
        print("overall model: ", str(self.model_overall), file=print_device)

    def plot_slice(self, metallicity_index, density_index, 
                   fitted=False, verbose=False):
        if verbose:
            print_device = sys.stdout
        else:
            print_device = open(os.devnull, mode="w")

        metallicity = self.metallicities_1D[metallicity_index]
        density     = self.densities_1D[density_index]
        mask = np.isclose(self.densities, density, atol=0) \
             & np.isclose(self.metallicities, metallicity, atol=0) \
             & self.usable & (self.momenta>0) # & ~self.excluded

        print(file=print_device)
        print("=================", file=print_device )
        print("metallicity: ", metallicity, file=print_device)
        print("density:     ", density, file=print_device)
        
        x_fit = np.logspace(0, 3.25, num=100)
        y_fit_overall = self.model_overall(metallicity, density, x_fit)

        if fitted:
            model_tmp = Momentum_Model(1e3, 1e4, 
                                       0, 0, 
                                       0, 0,
                                       2, -.12)

            fixed = np.array([False, False, 
                              True, True, 
                              True, True, 
                              False, True])
            
            x = (self.metallicities[mask], self.densities[mask], self.num_SNe[mask])
            y = self.momenta[mask] / (self.num_SNe[mask] * 100 * M_solar)

            fitted = False
            if len(x[0]) >= (~fixed).sum():
                fitted = True
                model_tmp.fit(x,y, fixed=fixed)
            models.append(model_tmp)
            fitted = False
            
            y_fit = model_tmp(metallicity, density, x_fit)
                                      
            print(str(models_tmp), file=print_device)

        plt.figure()
     
        print("number plotted: ", sum(mask), file=print_device)
        plt.scatter(self.num_SNe[mask], 
                    self.momenta[mask] / (self.num_SNe[mask] * 100 * M_solar),
                    marker= "o",
                    s=100,
                    label="data")
        
        plt.xlim(np.min(x_fit)/3 , np.max(x_fit)*3)
        
        plt.plot(x_fit, y_fit_overall, 
                 label="overall fit")
        if fitted:
            plt.plot(x_fit, y_fit,  label="fit", linestyle="--")
        plt.xscale("log")
        plt.legend(loc="best", frameon=True, shadow=True)
        plt.xlabel(r"$N_{SNe}$ ")
        plt.ylabel(r"$p / (100$ $M_\odot$ $N_\mathrm{SNe})$ $[\mathrm{km}$ $\mathrm{s}^{-1}]$")

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

    def print_latex_version(self, filename="plots_for_paper/equations/best_fit_equation.tex"):
        with open(filename, mode="w") as f:
            s = r"""\begin{{equation}}
        \frac{{p}}{{100 M_\odot N_\mathrm{{SNe}} }} = f\left( y_1, y_2 \right)
\end{{equation}}
        where
\begin{{align}}
        y_1 &= {0} \left(\frac{{Z}}{{Z_\odot}} \right)^{{{1:.2f}}} \left( \frac{{\rho}}{{m_p \text{{ cm}}^{{-3}}}}\right)^{{{2:.2f}}} \left( N_\mathrm{{SNe}} \right)^{{{3:.2f}}} \\
        y_2 &= {4} \left(\frac{{Z}}{{Z_\odot}} \right)^{{{5:.2f}}} \left( \frac{{\rho}}{{m_p \text{{ cm}}^{{-3}}}}\right)^{{{6:.2f}}} \left( \frac{{N_\mathrm{{SNe}}}}{{10^3}} \right)^{{{7:.2f}}}
\end{{align}}""".format(
                parse_into_scientific_notation(self.params[0]),
                self.params[2], 
                self.params[4], 
                self.params[6],
                parse_into_scientific_notation(self.params[1]), 
                self.params[3], 
                self.params[5], 
                self.params[7])
            print(s,file=f)



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
