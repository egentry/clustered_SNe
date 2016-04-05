from __future__ import print_function, division

from matplotlib import pyplot as plt


import numpy as np
import pandas as pd
from scipy import optimize, stats

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


class AggregatedResults(object):
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

    def get_MLE_fit(self, verbose=False):
        if verbose:
            print_device = sys.stdout
        else:
            print_device = open(os.devnull, mode="w")

        mask = self.usable & (self.momenta>0)
        xs = (self.metallicities[mask],
              self.densities[mask],
              self.num_SNe[mask])
        momenta = self.momenta[mask] / (self.num_SNe[mask] * 100 * M_solar)

        MLE_fit = MLEFit()
        MLE_fit.fit(self.metallicities[mask],
                    self.densities[mask],
                    self.num_SNe[mask],
                    momenta)

        return MLE_fit

    def get_Bayesian_fit(self, 
            generate_new_samples=False):

        Bayesian_fit = BayesianFit(self)
        if generate_new_samples:
            Bayesian_fit.generate_df()
        else:
            Bayesian_fit.load_df()

        return Bayesian_fit

    def plot_slice(self, metallicity, density,
                   with_MLE_fit = False, MLE_fit = None,
                   with_Bayesian_fit = False, Bayesian_fit = None,
                   verbose=False):
        """Plots a slice of momentum vs. N_SNe, at fixed density, metallicity

        Parameters
        ----------
        metallicity : float
            - gas metallicity (mass fraction)
        density : float
            - gas density [g cm**-3]
            
        with_MLE_fit : Optional[bool]
            - Should the plot include the MLE line-of-best fit?
            - If True, requires MLE_fit of type MLEFit
        MLE_fit : Optional[None or MLEFit (see with_MLE_fit)]

        with_Bayesian_fit : Optional[bool]
            - Should the plot include the Bayesian predictive median 
              and 1 sigma contours?
            - If True, requires Bayesian_fit of type BayesianFit
        Bayesian_fit : Optional[None or BayesianFit]

        verbose : Optional[bool]
            - print basic information about the plot
        """
        if verbose:
            print_device = sys.stdout
        else:
            print_device = open(os.devnull, mode="w")

        mask = np.isclose(self.densities, density, atol=0) \
             & np.isclose(self.metallicities, metallicity, atol=0) \
             & self.usable & (self.momenta>0) # & ~self.excluded

        print(file=print_device)
        print("=================",             file=print_device )
        print("metallicity:    ", metallicity, file=print_device)
        print("density:        ", density,     file=print_device)
        print("number plotted: ", sum(mask),   file=print_device)
        
        plt.figure()
     
        plt.scatter(self.num_SNe[mask], 
                    self.momenta[mask] / (self.num_SNe[mask] * 100 * M_solar),
                    marker= "o",
                    s=100,
                    label="data")

        if with_MLE_fit:
            x_fit = np.logspace(-.5, 3.5, num=100)
            y_fit = MLE_fit(metallicity, density, x_fit)
            plt.plot(x_fit, y_fit, label="MLE Fit")

        if with_Bayesian_fit:
            x_fit = np.logspace(-.5, 3.5, num=100)
            y_fit = Bayesian_fit.generate_predictive_percentiles(metallicity, density, x_fit)
            plt.plot(x_fit, y_fit[:,2], label="median( predictive )")
            plt.fill_between(x_fit, y_fit[:,1], y_fit[:,3], alpha=.5)
        
        plt.xlim(10**-.5 , 10**3.5)
        plt.ylim(ymin=0)
        
        plt.xscale("log")
        plt.legend(loc="best", frameon=True, shadow=True)
        plt.xlabel(r"$N_{SNe}$ ")
        plt.ylabel(r"$p / (100$ $M_\odot$ $N_\mathrm{SNe})$ $[\mathrm{km}$ $\mathrm{s}^{-1}]$")


#####################


def momentum_model(metallicity, density, N_SNe,
                   model_parameters):
    """A two-power law parameterization of our momentum results

    This function only *evaluates* the model.
    To fit the model, or save the best-fit parameters, see MLEFit and BayesianFit

    Parameters
    ----------
    metallicity : array-like, size N
        - gas metallicity (mass fraction)
    density : array-like, size N
        - gas density, in terms of [g cm**-3]
    N_SNe : array-like, size N
        - total number of SNe
    model_parameters : ModelParameters
        - ModelParameters.sigma_squared is not used here

    Returns
    -------
    momentum : array-like, size N
        - in the units of log10_norm_* (probably [km s**-1])
    """

    log10_norm_low,      log10_norm_high, \
    eta_metallicity_low, eta_metallicity_high, \
    eta_density_low,     eta_density_high, \
    eta_N_SNe_low,       eta_N_SNe_high = model_parameters.get_primary_parameters()
    
    norm_low  = 10**log10_norm_low
    norm_high = 10**log10_norm_high
            
    momentum_low  = norm_low  \
                * (metallicity / metallicity_solar)**eta_metallicity_low \
                * (density / m_proton)**eta_density_low \
                * (N_SNe / (1))**eta_N_SNe_low
    
    momentum_high = norm_high \
                * (metallicity / metallicity_solar)**eta_metallicity_high \
                * (density / m_proton)**eta_density_high \
                * (N_SNe / (1000))**eta_N_SNe_high
    
    momentum = (momentum_low * momentum_high) / (momentum_low + momentum_high)
        
    return momentum


#########################################


class ModelParameters(object):
    """Container for the model parameters of our 2 power law fit"""
    def __init__(self, 
                 log10_norm_low, log10_norm_high,
                 eta_metallicity_low, eta_metallicity_high,
                 eta_density_low, eta_density_high,
                 eta_N_SNe_low, eta_N_SNe_high,
                 sigma_squared=0):
        """
        Parameters
        (Note: *_low refers to low SNe power law, *_high refers to high SNe power law)
        All are type float
        ----------
        log10_norm_low  : normalization of power law
        log10_norm_high : normalization of power law
        eta_metallicity_low  : power law index on metallicity dependence
        eta_metallicity_high : power law index on metallicity dependence
        eta_density_low  : power law index on ISM density dependence
        eta_density_high : power law index on ISM density dependence
        eta_N_SNe_low  : power law index on N_SNe dependence
        eta_N_SNe_high : power law index on N_SNe dependence
        
        sigma_squared : [Optional] variance of uncertainty added to model"""
                
        self.log10_norm_low  = log10_norm_low
        self.log10_norm_high = log10_norm_high
        self.eta_metallicity_low  = eta_metallicity_low
        self.eta_metallicity_high = eta_metallicity_high
        self.eta_density_low  = eta_density_low
        self.eta_density_high = eta_density_high
        self.eta_N_SNe_low  = eta_N_SNe_low
        self.eta_N_SNe_high = eta_N_SNe_high
        self.sigma_squared = sigma_squared
        
    def get_primary_parameters(self):
        return [self.log10_norm_low,      self.log10_norm_high,
                self.eta_metallicity_low, self.eta_metallicity_high,
                self.eta_density_low,     self.eta_density_high,
                self.eta_N_SNe_low,       self.eta_N_SNe_high]

    def get_nuisance_parameters(self):
        return [self.sigma_squared,]

    def get_parameters(self):
        return [*self.get_primary_parameters(), *self.get_nuisance_parameters()]

    def copy(self):
        return self.__class__(*self.get_parameters())


#########################################


class MLEFit(object):
    """Maximum Likelihood Estimated (MLE) fit to our momentum_model"""

    def __init__(self,
                model_parameters=ModelParameters(
                    np.log10(1e3), np.log10(1e4), # log10_norm
                    0,             0,             # eta_metallicity
                    .25,           .2,            # eta_density
                    2,             -.08           # eta_N
                )):
        """
        Parameters
        (Note: *_low refers to low SNe power law, *_high refers to high SNe power law)
        All are type float; All are initial guesses
        ----------
        model_parameters : ModelParameters
            - ModelParameters.sigma_squared is not used here
        """
        
        params_0 = model_parameters
        
        self.model_parameters_initial = model_parameters
        self.model_parameters         = model_parameters.copy()
    
    def __call__(self, metallicity, density, N_SNe):        
        return momentum_model(metallicity, density, N_SNe,
                              self.model_parameters)
        
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
        s = s.format(10**self.model_parameters.log10_norm_low, 
                         self.model_parameters.eta_metallicity_low, 
                         self.model_parameters.eta_density_low, 
                         self.model_parameters.eta_N_SNe_low,
                     10**self.model_parameters.log10_norm_high, 
                         self.model_parameters.eta_metallicity_high, 
                         self.model_parameters.eta_density_high, 
                         self.model_parameters.eta_N_SNe_high)
        
        return s

    def print_latex_version(self, filename="plots_for_paper/equations/best_fit_equation.tex"):
        with open(filename, mode="w") as f:
            s = r"""\begin{{equation}}
        \frac{{p}}{{100 M_\odot N_\mathrm{{SNe}} }} = f\left( y_1, y_2 \right)
\end{{equation}}
        where
\begin{{align}}
        y_1 &= {0} \left(\frac{{Z}}{{Z_\odot}} \right)^{{{1:.2f}}} \left( \frac{{\rho}}{{m_p \text{{ cm}}^{{-3}}}}\right)^{{{2:.2f}}} \left( N_\mathrm{{SNe}} \right)^{{{3:.2f}}} \label{{eq:model:fewSNe}} \\
        y_2 &= {4} \left(\frac{{Z}}{{Z_\odot}} \right)^{{{5:.2f}}} \left( \frac{{\rho}}{{m_p \text{{ cm}}^{{-3}}}}\right)^{{{6:.2f}}} \left( \frac{{N_\mathrm{{SNe}}}}{{10^3}} \right)^{{{7:.2f}}} \label{{eq:model:manySNe}}
\end{{align}}""".format(
                parse_into_scientific_notation(10**self.model_parameters.log10_norm_low),
                self.model_parameters.eta_metallicity_low, 
                self.model_parameters.eta_density_low, 
                self.model_parameters.eta_N_SNe_low,
                parse_into_scientific_notation(10**self.model_parameters.log10_norm_high), 
                self.model_parameters.eta_metallicity_high, 
                self.model_parameters.eta_density_high, 
                self.model_parameters.eta_N_SNe_high)
            print(s, file=f)



    def fit(self, metallicity, density, N_SNe, momenta, fixed=None, **kwargs):
        """Fits the momentum using the initial self.model_parameters_initial;
        Saves the best-fit parameters into self.model_parameters

        Parameters
        ----------
            metallicity : array-like, size N
                - gas metallicity (mass fraction)
            density : array-like, size N
                - gas density, in terms of [g cm**-3]
            N_SNe : array-like, size N
                - total number of SNe
            momenta : momentum efficiency (observed), array-like, size N
                - momentum / (100*M_solar*N_SNe) in units of [km s**-1]
            fixed : mask of values to be held fixed; None or np.ndarray (dtype=bool)
                - fixed values set to True
                - free  values set to False
                - None for all values free
            **kwargs : to be passed to optimize.curve_fit; dict
                - does NOT get passed to the fitted function

        Side effects
        ------------
        Overwrites self.model_parameters with latest fit results
        
        Notes
        -----
        Fits momentum_model, using a simple least-squares frequentist method
        """
        p0 = np.array(self.model_parameters_initial.get_primary_parameters())
        
        if fixed is None:
            fixed = np.full_like(p0, False, dtype=bool)
        elif isinstance(fixed, np.ndarray):
            if not fixed.dtype == bool:
                raise TypeError("fixed must be an np.ndarray of bools")
        else:
            raise TypeError("fixed must be type None or np.ndarray")
            
        def _f(x, *args):
            metallicity, density, N_SNe = x
            p0[~fixed] = args
            return momentum_model(metallicity, density, N_SNe,
                                  ModelParameters(*p0))
        
        popt, _ = optimize.curve_fit(_f, (metallicity, density, N_SNe),
                                     momenta, 
                                     p0=p0[~fixed],
                                     xtol=.01,
                                     **kwargs)
        p0[~fixed] = popt
        self.model_parameters = ModelParameters(*p0)


#########################################


class BayesianFit(object):
    """Creates, and holds the Bayesian fit to our momentum_model"""
    def __init__(self, aggregated_results):
        self.aggregated_results = aggregated_results
        
        self.mask =  self.aggregated_results.usable & (self.aggregated_results.momenta>0)

        self.model_inputs = (self.aggregated_results.metallicities[self.mask],
                             self.aggregated_results.densities[self.mask],
                             self.aggregated_results.num_SNe[self.mask])
        
        self.momentum_measured = self.aggregated_results.momenta[self.mask] \
            / (self.aggregated_results.num_SNe[self.mask] * 100 * M_solar)
        
        self.n_measurements = self.momentum_measured.size 
        
        self.generate_initial_parameters()
        
    def save_df(self, filename="posterior_samples.h5"):
        with pd.HDFStore(filename) as store:
            store["df"] = self.df
    
    def load_df(self, filename="posterior_samples.h5"):
        with pd.HDFStore(filename) as store:
            self.df = store["df"]
    
    def generate_initial_parameters(self):
        """Use the MLE parameters, and estimate sigma_squared from the variance of the errors"""
        MLE_fit = self.aggregated_results.get_MLE_fit()
        self.model_parameters_initial = MLE_fit.model_parameters

        momentum_modeled = MLE_fit(*self.model_inputs)
        sigma_squared = np.var(momentum_modeled - self.momentum_measured)

        self.model_parameters_initial.sigma_squared = sigma_squared 
        
    def log_prior(self, model_parameters):
        return -np.log(model_parameters.sigma_squared)
    
    def log_likelihood(self, model_parameters):
        momentum_modeled  = momentum_model(*self.model_inputs,
                                           model_parameters)

        return np.sum(- .5 * (momentum_modeled - self.momentum_measured)**2\
                              /model_parameters.sigma_squared \
                      - .5 * np.log(model_parameters.sigma_squared))
    
    def log_posterior(self, model_parameters):
        return self.log_prior(model_parameters) + self.log_likelihood(model_parameters)
    
    def draw_new_samples(self, model_parameters):

        # Use Metropolis-Hastings for primary parameters
        step_size = [.05, .05, # log10_norm_*
                     .05, .05, # eta_metallicity_*
                     .05, .05, # eta_density_*
                     .5, .05,  # eta_N_SNe_*
                     0]

        acceptances = np.empty_like(model_parameters.get_parameters(), dtype=int)
        model_parameters_new = model_parameters.copy()

        for i, parameter_i in enumerate(model_parameters.get_parameters()[:-1]):
            parameters_tmp = np.array(model_parameters_new.get_parameters())
            parameters_tmp[i] = stats.norm.rvs(loc=parameter_i, scale=step_size[i])
            model_parameters_tmp = ModelParameters(*parameters_tmp)

            log_prob_ratio =    self.log_posterior(model_parameters_tmp) \
                              - self.log_posterior(model_parameters_new)

            if log_prob_ratio > np.log(stats.uniform.rvs()):
                model_parameters_new = model_parameters_tmp
                acceptances[i] = 1
            else:
                acceptances[i] = 0
        
        
        # Gibbs sample a new sigma_squared value
        momentum_modeled  = momentum_model(*self.model_inputs,
                                           model_parameters_new)
        shape = self.n_measurements/2
        scale = np.sum(.5*(momentum_modeled - self.momentum_measured)**2)
        sigma_squared_new = stats.invgamma.rvs(shape, scale=scale)
        model_parameters_new = ModelParameters(*model_parameters_new.get_primary_parameters(),
                                               sigma_squared_new)
        acceptances[-1] = 1

        return model_parameters_new, acceptances

    def generate_df(self, n_samples=100000, verbose=False):
        n_burn  = n_samples//10
        n_draws = n_samples + n_burn
        
        samples = np.empty((n_draws, 
                            len(self.model_parameters_initial.get_parameters())))
        samples[0] = self.model_parameters_initial.get_parameters()
        
        acceptances_total = np.zeros(len(self.model_parameters_initial.get_parameters()),
                                     dtype=int)
        for i in range(1, n_draws):
            model_parameters_new, acceptances = self.draw_new_samples(ModelParameters(*samples[i-1]))
            
            samples[i] = model_parameters_new.get_parameters()
            
            if i >= n_burn:
                acceptances_total += acceptances
        
        samples = samples[n_burn:]
        
        columns = ["log10_norm_low",      "log10_norm_high",
                   "eta_metallicity_low", "eta_metallicity_high",
                   "eta_density_low",     "eta_density_high",
                   "eta_N_SNe_low",       "eta_N_SNe_high",
                   "sigma_squared"]
        
        if verbose:
            print(acceptances_total/n_samples)
                
        self.df =  pd.DataFrame(samples, columns=columns)
        
    def create_trace_plots(self):
        for i, column in enumerate(self.df):
            plt.figure()
            self.df[column].plot(label="MCMC")
            plt.ylabel(column)
            plt.title("Trace Plot")
            plt.axhline(self.model_parameters_initial.get_parameters()[i], c="r", label="MLE")
            plt.legend(loc="best")
    
    def create_autocorrelation_plots(self):
        for column in self.df:
            plt.figure()
            plt.acorr(self.df[column], 
                      maxlags=min(1000, self.df.shape[0]), usevlines=False, linestyle="-", label=column)
            plt.legend(loc="best")
            plt.xlabel("Lag")
            plt.ylabel("Autocorrelation")
            
    def create_corner_plots(self):
        import corner
        labels = [r"$\log p_\mathrm{low}$", r"$\log p_\mathrm{high}$",
          r"$\eta_{Z,low}$", r"$\eta_{Z, high}$",
          r"$\eta_{\rho, low}$", r"$\eta_{\rho, high}$",
          r"$\eta_{N,low}$", r"$\eta_{N, high}$",
          r"$\sigma^2$"]
        figure = corner.corner(self.df, 
                       label_kwargs={"fontsize":40},
                       labels=labels,
                       show_titles = True)
        return figure
    
    def generate_predictive_percentiles(self, metallicity, density, num_SNe, 
                                        percentiles=np.array([2.28, 15.87, 50, 84.13, 97.72]),
                                        subsample_factor = 100, num_noise_realizations=100):
        """For the given metallicity, density, SNe this generates the predictive distribution,
        and returns the momentum corresponding to a certain percentile
        
        Parameters
        ----------
        metallicity : 1D array_like, size N or 1 
            - gas mass fraction
        density : 1D array_like, size N or 1 
            - units of [g cm**-3]
        num_SNe : 1D array_like, size N or 1 
        percentiles : Optional[array-like (Mx1)]
            - return the momenta corresponding to these percentiles for each metallicity,density,num_SNe
            - default relative to median: -2 sigma, -1 sigma, 0, +1 sigma, +2 sigma
        subsample_factor : Optional[int]
            - when iterating through the MCMC samples, take steps of size subsample_factor
        num_noise_realizations : Optional[int]
            - for a given MCMC sample, we need to create `num_noise_realizations` of the noise
        
        Returns
        -------
        momenta : arraylike (NxM)
            - momentum / (100 * M_solar * num_SNe)
            - units of [km s**-1] momentum / (100 * M_solar * num_SNe)
        """
        metallicity = np.array(metallicity).flatten()
        density     = np.array(density).flatten()
        num_SNe     = np.array(num_SNe).flatten()
        
        N = max(metallicity.size, density.size, num_SNe.size)
        
        metallicity = metallicity * np.ones(N)
        density     = density     * np.ones(N)
        num_SNe     = num_SNe     * np.ones(N)
            
        M = len(percentiles)
        
        df_tmp = self.df.loc[::subsample_factor]
        J = df_tmp.shape[0]
        
        momenta = np.empty((N,M))
        for i in range(N):            
            ys = np.empty((J, num_noise_realizations))
            for j, row in enumerate(df_tmp.iterrows()):
                parameters = list(row[1])
                sigma_squared = parameters[8]
                ys[j] = momentum_model(metallicity[i], density[i], num_SNe[i],
                                       ModelParameters(*parameters))
                ys[j] += stats.norm.rvs(scale=np.sqrt(sigma_squared), size=num_noise_realizations)
            momenta[i] = np.percentile(ys, percentiles)
                
        return momenta
        

        