
import numpy as np
import pandas as pd

import glob

## Boilerplate path hack to give access to full clustered_SNe package
import sys, os
if __package__ is None:
    if os.pardir not in sys.path[0]:
        file_dir = os.getcwd()
        sys.path.insert(0, os.path.join(file_dir, 
                                        os.pardir, 
                                        os.pardir))

from clustered_SNe.analysis.constants import m_proton, metallicity_solar
from clustered_SNe.analysis.parse import RunSummary, Overview



from astropy.modeling import models, fitting


from IPython.display import display, Math

import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext, \
                              LogLocator

energies = ["tot", "kin", "int" ]

#####################


class ThorntonParameterStudy(object):
    """
    Holds the reduced and fitted results for a full parameter study 
    in the style of Thornton et al. (1998)

    """

    formatter = LogFormatterMathtext(labelOnlyBase=False)
    locator   = LogLocator(subs=np.logspace(0,1,num=5,endpoint=False)[1:])

    metallicity_fit_label    = r"Fit($Z$)"
    number_density_fit_label = r"Fit($n$)"
    simultaneous_fit_label   = r"Fit($Z$, $n$)"
    my_label = "Numeric Results"
    Thornton_label = "Thornton et al. (1998)"



    def __init__(self, data_dir,
                       with_cooling = True,
                       verbose = False):
        self.data_dir = data_dir
        
        columns = ["background_density", "background_temperature", "metallicity",
                   "number_density",
                   "E_R_kin", "E_R_int", "E_R_tot", "t_f", "R_shock"]

        overview_filenames = glob.glob(os.path.join(data_dir, "*_overview.dat"))

        overviews = [Overview(filename) for filename in overview_filenames]
        overviews = [overview for overview in overviews if overview.with_cooling == with_cooling]

        self.df = pd.DataFrame(index=range(len(overviews)), columns=columns, dtype=float)

        for i, overview in enumerate(overviews):
            
            run_summary = RunSummary(data_dir=data_dir, id=overview.id)
            cooling_index = np.argmin(np.abs(run_summary.times - run_summary.t_f))

            self.df.loc[i] = [overview.background_density, 
                         overview.background_temperature,
                         overview.metallicity,
                         overview.background_density / m_proton,
                         run_summary.E_R_kin[cooling_index], 
                         run_summary.E_R_int[cooling_index], 
                         run_summary.E_R_tot[cooling_index], 
                         run_summary.t_f,
                         run_summary.R_shock[cooling_index]]

            if verbose:
                print()
                print("metallicity: ", overview.metallicity)
                print("background density:", overview.background_density)
                print("cooling index: ", cooling_index)



        ## add models
        sorted_args = np.argsort(self.df.number_density)
        self.slice_by_metallicity = sorted_args[np.isclose(self.df.metallicity[sorted_args], metallicity_solar, atol=0)]

        sorted_args = np.argsort(self.df.metallicity)
        self.slice_by_number_density = sorted_args[np.isclose(self.df.number_density[sorted_args], .1 * 1.33, atol=0)]

        self.thornton_models       = dict()

        self.simultaneous_models   = dict()
        self.number_density_models = dict()
        self.metallicity_models    = dict()

        for energy in energies:
            self.thornton_models[energy] = PublishedModel(energy)
            self.simultaneous_models[energy] = SimultaneousModel(self.df.number_density, 
                                                              self.df.metallicity, 
                                                              self.df["E_R_{0}".format(energy)])
            self.number_density_models[energy] = NumberDensityModel(self.df.number_density[self.slice_by_metallicity],
                                                               self.df["E_R_{0}".format(energy)][self.slice_by_metallicity])
            self.metallicity_models[energy] = MetallicityModel(self.df.metallicity[self.slice_by_number_density],
                                                          self.df["E_R_{0}".format(energy)][self.slice_by_number_density])

        if verbose:
            for energy in energies:

                print("")
                print("Thornton Model, E_R_{0}:".format(energy))
                self.thornton_models[energy].print_model()


            print("========================")

            for energy in energies:
                print("")
                print("simultaneous model, E_R_{0}".format(energy))
                self.simultaneous_models[energy].print_model()




            print("========================")

            for energy in energies:
                print("")
                print("number_density_fit_E_R_{0}: ".format(energy))
                self.number_density_models[energy].print_model()


            print("========================")

            for energy in energies:
                print("")
                print("metallicity_fit_E_R_{0}: ".format(energy))
                self.metallicity_models[energy].print_model()

    def plot_one_number_density(self, energy="tot"):
        mask = self.slice_by_metallicity

        y_variable = "E_R_" + energy

        self.number_density_models[energy].print_model()

        plt.figure()
        plt.plot(self.df.number_density[mask], self.df[y_variable][mask],
                 marker="o", linestyle="",
                 label=self.my_label)
        plt.plot(self.df.number_density[mask],
                 self.number_density_models[energy](self.df.number_density[mask]),
                 label=self.number_density_fit_label)
        plt.plot(self.df.number_density[mask],
                 self.simultaneous_models[energy](self.df.metallicity[mask],
                                             self.df.number_density[mask]),
                 label=self.simultaneous_fit_label)
        plt.plot(self.df.number_density[mask],
                 self.thornton_models[energy](self.df.metallicity[mask],
                                self.df.number_density[mask]),
                 label=self.Thornton_label, linestyle="--")
        plt.legend(loc="best")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel(r"$n$ [cm$^{-3}$]")
        plt.ylabel(r"$E_{{R, \mathrm{{ {0} }} }}$ [ergs]".format(energy))

        xmin, xmax = plt.xlim()
        plt.xlim(xmin/5, xmax*5)

        ax = plt.gca()
        ax.yaxis.set_minor_formatter(self.formatter)
        ax.yaxis.set_minor_locator(self.locator)


    def plot_one_metallicity(self, energy="tot"):
        mask = self.slice_by_number_density

        y_variable = "E_R_" + energy
        
        self.metallicity_models[energy].print_model()
        
        plt.figure()
        plt.plot(self.df.metallicity[mask] / metallicity_solar,
                 self.df[y_variable][mask],
                 marker="o", linestyle="", label="Numeric Results")
        plt.plot(self.df.metallicity[mask] / metallicity_solar,
                 self.metallicity_models[energy](self.df.metallicity[mask]),
                 label=self.metallicity_fit_label)
        plt.plot(self.df.metallicity[mask] / metallicity_solar,
                 self.simultaneous_models[energy](self.df.metallicity[mask],
                                             self.df.number_density[mask]),
                 label=self.simultaneous_fit_label)
        plt.plot(self.df.metallicity[mask] / metallicity_solar,
                 self.thornton_models[energy](self.df.metallicity[mask],
                                self.df.number_density[mask]),
                 label=self.Thornton_label, linestyle="--")
        plt.legend(loc="best")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel(r"$\log[Z / Z_\odot]$")
        plt.ylabel(r"$E_{{R, {0} }}$ [ergs]".format(energy))

        xmin, xmax = plt.xlim()
        plt.xlim(xmin/5, xmax*5)
        
        ax = plt.gca()
        ax.yaxis.set_minor_formatter(self.formatter)
        ax.yaxis.set_minor_locator(self.locator)


#####################


# astropy's models don't play well with fitting power laws
# at least as far as making ones that look good in loglog plots
# 
# So instead, I'm going to go with these ugly wrapper classes

class SimultaneousModel(object):
    """Fit both metallicity and number density dependence"""
    def __init__(self, metallicities, number_densities, energies):
        super(SimultaneousModel, self).__init__()
        self.model = self.fit(metallicities, number_densities, energies)
        # should I keep the input metallicities, etc as attributes?
        
        
    def __call__(self, metallicities, number_densities):
        log_number_densities = np.log10(number_densities)
        log_metallicities    = np.log10(metallicities / metallicity_solar)
        
        log_energies = self.model(log_metallicities, log_number_densities)
        
        return 10**log_energies
        
    def fit(self, metallicities, number_densities, energies):
        fitter = fitting.LinearLSQFitter()
        log_number_densities = np.log10(number_densities)
        log_metallicities    = np.log10(metallicities / metallicity_solar)
        log_energies         = np.log10(energies)
        
        model = models.Polynomial2D(1)
        
        fitted = fitter(model, log_metallicities,
                               log_number_densities,
                               log_energies)
        return fitted
    
    def print_model(self):
        result  = r"E = "
        result += r"{0:.2f} \cdot 10^{{49}}".format(10**self.model.c0_0.value/1e49)
        result += r" \cdot (Z / Z_\odot)^{{ {0:.2f} }}".format(self.model.c1_0.value)
        result += r" \cdot n_0^{{ {0:.2f} }}".format(self.model.c0_1.value)
        result += r" \text{ [ergs]}"
        display(Math(result))    
        
        

class NumberDensityModel(object):
    """Eventually I should make this a real astropy model subclass"""
    def __init__(self, number_densities, energies):
        super(NumberDensityModel, self).__init__()
        self.model  = self.fit(number_densities, energies)
        
    def __call__(self, number_densities):
        log_number_densities = np.log10(number_densities)

        log_energies = self.model(log_number_densities)
        return 10**log_energies

    def fit(self, number_densities, energies):
        fitter = fitting.LinearLSQFitter()
        log_number_densities = np.log10(number_densities)
        log_energies         = np.log10(energies)

        log_number_density_model = models.Linear1D()
        fitted = fitter(log_number_density_model,
                        log_number_densities,
                        log_energies)
        return fitted

    def print_model(self):
        result  = r"E = "
        result += r"{0:.2f} \cdot 10^{{49}}".format(10**self.model.intercept.value/1e49)
        result += r" \cdot n_0^{{ {0:.2f} }}".format(self.model.slope.value)
        result += r" \text{ [ergs]}"
        display(Math(result))    
    
    
    
class MetallicityModel(object):
    """Eventually I should make this a real astropy model subclass"""
    def __init__(self, metallicities, energies):
        super(MetallicityModel, self).__init__()
        self.model  = self.fit(metallicities, energies)
        
    def __call__(self, metallicities):
        log_metallicities = np.log10(metallicities / metallicity_solar)

        log_energies = self.model(log_metallicities)
        return 10**log_energies

    def fit(self, metallicities, energies):
        fitter = fitting.LinearLSQFitter()
        log_metallicities = np.log10(metallicities / metallicity_solar)
        log_energies      = np.log10(energies)

        log_metallicity_model = models.Linear1D()
        fitted = fitter(log_metallicity_model,
                        log_metallicities,
                        log_energies)
        return fitted

    def print_model(self):
        result  = r"E = "
        result += r"{0:.2f} \cdot 10^{{49}}".format(10**self.model.intercept.value/1e49)
        result += r" \cdot (Z / Z_\odot)^{{ {0:.2f} }}".format(self.model.slope.value)
        result += r" \text{ [ergs]}"
        display(Math(result))
        


class PublishedModel(object):
    def __init__(self, energy):
        super(PublishedModel, self).__init__()
        
        if energy in energies:
            self.energy = energy
        else:
            raise ValueError("Incorrect input of 'energy': " + energy)
        
    def __call__(self, metallicity, number_density):
        if self.energy is "tot":
            return self.E_tot(metallicity, number_density)
        elif self.energy is "kin":
            return self.E_kin(metallicity, number_density)
        elif self.energy is "int":
            return self.E_int(metallicity, number_density)
    
    def E_tot(self, metallicity, number_density):
        E_tot  = self.E_kin(metallicity, number_density)
        E_tot += self.E_int(metallicity, number_density)
        return E_tot

    def E_kin(self, metallicity, number_density):
        # the awkward powers of 0 help us return
        # an array of the correct shape
        E_kin = 8.52e49 * metallicity**0 * number_density**0
        return E_kin
    
    def E_int(self, metallicity, number_density):
        E_int = 1.83e49 * number_density**-.23 \
                * (metallicity / metallicity_solar)**-.24
        return E_int
    
    def print_model(self):
        if self.energy is "tot":
            result = r"E_{tot} = E_{int} + E_{kin}"
        elif self.energy is "kin":
            result  = r"E_{kin} = 8.52 \cdot 10^{49}"
            result += r" \text{ [ergs]}"
        elif self.energy is "int":
            result  = r"E_{int} = 1.83 \cdot 10^{49}"
            result += r" \cdot (Z / Z_\odot)^{-0.24}"
            result += r" n_0^{-0.23}" 
            result += r" \text{ [ergs]}"

        display(Math(result))   


