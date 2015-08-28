import os
import glob

import numpy as np
import pandas as pd

from astropy.convolution import convolve, Gaussian1DKernel

## Boilerplate path hack to give access to full SNe package
import sys
if __package__ is None:
    if os.pardir not in sys.path[0]:
        file_dir = os.path.dirname(__file__)
        sys.path.insert(0, os.path.join(file_dir, 
                                        os.pardir, 
                                        os.pardir))

from SNe.analysis.constants import m_proton, pc, M_solar, gamma, E_0, \
                                   metallicity_solar, yr

from SNe.analysis.helper_functions import calculate_mean_molecular_weight, \
                                          calculate_mass, \
                                          calculate_kinetic_energy, \
                                          calculate_internal_energy, \
                                          calculate_momentum, \
                                          calculate_c_ad, \
                                          calculate_entropy, \
                                          calculate_temperature, \
                                          calculate_w_cell, \
                                          calculate_crossing_time


class RunSummary(object):
    """ 
        Allows return variables to be stored in function arguments
        
        The better way to do this would be to get the return variables
        from the functions I pass to interact*
        
        But I haven't figured out how to do that yet.

        Other attributes will be added
    """
    def __init__(self, df=None, zones=None, E_tot=None, E_int=None, E_kin=None):
        super(RunSummary, self).__init__()
        self.df    = df
        self.zones = zones
        self.E_tot = E_tot
        self.E_int = E_int
        self.E_kin = E_kin



class Overview(object):
    """Basic overview of a given ./SNe run


    Attributes
    ----------
    id : str 
    metallicity : Optional[float]
    background_density : Optional[float]
    background_temperature : Optional[float]
    with_cooling : Optional[bool]



    """
    def __init__(self, filename):
        """Create an Overview object using an "overview.dat" style filename (see Output/ascii.c)

        Parameters
        ----------
        filename - should be a valid "overview.dat" style filename (see Output/ascii.c)

        Notes
        -----
        We can currently parse the following attributes:
            Metallicity : float
                [ mass fraction ]
            Background Density : float
                [ g cm^-3 ]
            Background Temperature : float
                [ K ]
            With Cooling : bool
            Number of SNe : int
            Cluster Mass : float
                [g]

        """
        super(Overview, self).__init__()
        
        self.id = os.path.basename(filename).split("_")[0]  
        self.dirname = os.path.dirname(filename)
        # Add trailing slash (if dirname isn't empty)
        self.dirname = os.path.join(self.dirname, "")
        
         # default, since earlier runs won't have this saved
        self.num_SNe = 0
        self.cluster_mass = 0


        f = open(filename, "r")
        for line in f:
            if "Metallicity" in line:
                self.metallicity = float(line.split()[1])
            elif "Background Density" in line:
                self.background_density = float(line.split()[2])
            elif "Background Temperature" in line:
                self.background_temperature = float(line.split()[2])
            elif "With cooling" in line:
                self.with_cooling = bool(int(line.split()[2]))
            elif "Number of SNe" in line:
                self.num_SNe = int(line.split()[-1])
            elif "Cluster Mass" in line:
                self.cluster_mass = float(line.split()[-1]) * M_solar
        f.close()

        SNe_filename = filename.replace("overview", "SNe") 
        if os.path.exists(SNe_filename):
            SNe = np.loadtxt(SNe_filename)
            if (SNe.shape[0] != self.num_SNe):
                raise ValueError("Number of SNe in datafile " +
                                 "doesn't match number listed in overview file")
            self.SNe_times          = SNe[:,0]
            self.SNe_initial_mass   = SNe[:,1]
            self.SNe_ejecta_mass    = SNe[:,2]
            self.SNe_ejecta_mass_Z  = SNe[:,3]
        else:
            self.SNe_times          = np.array([0.])
            self.SNe_initial_mass   = np.array([0.])
            self.SNe_ejecta_mass    = np.array([0.])
            self.SNe_ejecta_mass_Z  = np.array([0.])

        sorted_indices = np.argsort(self.SNe_times)
        self.SNe_times          = self.SNe_times[         sorted_indices]
        self.SNe_initial_mass   = self.SNe_initial_mass[  sorted_indices]
        self.SNe_ejecta_mass    = self.SNe_ejecta_mass[   sorted_indices]
        self.SNe_ejecta_mass_Z  = self.SNe_ejecta_mass_Z[ sorted_indices]

        return

    def __str__(self):
        string  = "id \t\t\t = {0}".format(self.id) + "\n"
        string += "metallicity \t\t = {0} ".format(self.metallicity) + "\n"
        string += "background density \t = {0} [g cm^-3]".format(self.background_density) + "\n"
        string += "background temperature \t = {0:e} [K]".format(self.background_temperature)
        return string


class ParseResults(object):
    """ 
        This should be taken out when I have more time
    """
    def __init__(self, checkpoint_filenames):
        super(ParseResults, self).__init__()
        self.checkpoint_filenames   = checkpoint_filenames 
        return


cols = ["Radius", "dR", "dV", "Density", 
        "Pressure", "Velocity", "Z", 
        "Temperature", "Energy", "Entropy", 
        "Mass", "M_int", "C_ad", "Crossing_time"]
cols_in   = cols[:-7]




def parse_run(data_dir="", id="", last_run=None):
    #this whole thing is a mess, and needs to be refactored

    checkpoint_filenames = glob.glob(os.path.join(data_dir,id + "*checkpoint_*.dat"))
    checkpoint_filenames = sorted(checkpoint_filenames)
    num_checkpoints = len(checkpoint_filenames)

    if num_checkpoints == 0:
        raise FileNotFoundError("No checkpoints found")
    
    # ensure that the id is actually the FULL id
    basename = os.path.basename(checkpoint_filenames[0])
    id = basename.split("checkpoint")[0]

    times = np.empty(num_checkpoints)
    for k, checkpoint_filename in enumerate(checkpoint_filenames):
        f = open(checkpoint_filename, 'r')
        line = f.readline()
        times[k] = float(line.split()[3])
        f.close()

    overview = Overview(os.path.join(data_dir, id + "overview.dat"))

    mu = calculate_mean_molecular_weight(overview.metallicity)

    E_int    = np.empty(num_checkpoints)
    E_kin    = np.empty(num_checkpoints)
    momentum = np.empty(num_checkpoints)
    E_tot    = np.empty(num_checkpoints)
    
    E_R_int  = np.empty(num_checkpoints) # Energy of the remnant
    E_R_kin  = np.empty(num_checkpoints) # Energy of the remnant
    E_R_tot  = np.empty(num_checkpoints) # Energy of the remnant
    R_shock  = np.empty(num_checkpoints) # size of the shock/remnant
    
    M_tot    = np.empty(num_checkpoints)
    Z_tot    = np.empty(num_checkpoints)
    zones    = np.empty(num_checkpoints)
    
    Luminosity = np.empty(num_checkpoints)


    #### PARSE DATAFILES INTO DATAFRAME
    df = pd.DataFrame()
    for k, checkpoint_filename in enumerate(checkpoint_filenames):
        array_tmp = np.loadtxt(checkpoint_filename, usecols=range(len(cols_in)))
        array_tmp = array_tmp[1:-1] # ignore guard cells
        index     = pd.MultiIndex.from_product([k, np.arange(array_tmp.shape[0])],
                                               names=["k","i"])
        df_tmp    = pd.DataFrame(array_tmp, columns=cols_in, index = index)

        df_tmp["Mass"]        = calculate_mass(df_tmp.Density.values,
                                               df_tmp.dV.values)
        df_tmp["M_int"]       = df_tmp.Mass.cumsum()
        df_tmp["Temperature"] = calculate_temperature(df_tmp.Pressure.values,
                                                      df_tmp.Density.values, mu)
        df_tmp["Energy"]      = calculate_internal_energy(df_tmp.Mass.values,
                                                          df_tmp.Pressure.values,
                                                          df_tmp.Density.values) \
                                                           / df_tmp.Mass.values
        df_tmp["Entropy"]     = calculate_entropy(df_tmp.Temperature.values, 
                                                  df_tmp.Density.values, mu)
        df_tmp["C_ad"]        = calculate_c_ad(df_tmp.Pressure.values,
                                               df_tmp.Density.values)
        w_cell = calculate_w_cell(df_tmp.Velocity.values)
        df_tmp["Crossing_time"] = calculate_crossing_time(df_tmp.C_ad.values,
                                                          df_tmp.Velocity.values,
                                                          w_cell,
                                                          df_tmp.dR.values )
        df_tmp["zones"]       = np.arange(array_tmp.shape[0])
        
        E_kin[k]    = calculate_kinetic_energy(df_tmp.Mass.values,
                                               df_tmp.Velocity.values)
        E_int[k]    = calculate_internal_energy(df_tmp.Mass.values,
                                                df_tmp.Pressure.values,
                                                df_tmp.Density.values)
        momentum[k] = calculate_momentum(df_tmp.Mass.values,
                                         df_tmp.Velocity.values)
        E_tot[k]    = E_kin[k] + E_int[k]
        M_tot[k]    = df_tmp.M_int[-1]
        Z_tot[k]    = np.sum(df_tmp.Mass * df_tmp.Z)
        zones[k]    = df_tmp.shape[0]
        
        # over_dense = df_tmp.Density != background_density
        over_dense = df_tmp.Density > overview.background_density * 1.0001
        if np.any(over_dense):
            R_shock[k] = np.max(df_tmp.Radius[over_dense])
        else:    
            R_shock[k] = df_tmp.Radius.iloc[0]
        remnant = df_tmp.Radius <= R_shock[k]
        E_R_kin[k] = calculate_kinetic_energy(df_tmp.Mass[remnant].values,
                                              df_tmp.Velocity[remnant].values)
        E_R_int[k] = calculate_internal_energy(df_tmp.Mass[remnant].values, 
                                               df_tmp.Pressure[remnant].values,
                                               df_tmp.Density[remnant].values)
        E_R_tot[k] = E_R_int[k] + E_R_kin[k]
        
        if k is 0:
            Luminosity[k] = 0
        else:
            Luminosity[k] = -(E_R_tot[k] - E_R_tot[k-1]) / (times[k] - times[k-1])

        df = pd.concat([df, df_tmp])
        
    df.Radius /= pc
    df.Mass   /= M_solar
    df.M_int  /= M_solar
    

    if type(last_run) is not RunSummary:
        # only overwrite last_run if it's a RunSummary
        # this logic is sketchy
        last_run= RunSummary()

    last_run.df         = df
    last_run.times      = times
    last_run.zones      = zones
    last_run.E_tot      = E_tot
    last_run.Z_tot      = Z_tot
    last_run.E_int      = E_int
    last_run.E_kin      = E_kin
    last_run.E_R_int    = E_R_int
    last_run.E_R_kin    = E_R_kin
    last_run.E_R_tot    = E_R_tot
    last_run.R_shock    = R_shock
    last_run.momentum   = momentum
    last_run.Luminosity = Luminosity

    last_run.overview = overview
    
    # filter for when initial transients have settled
    # assume that the settling time scales with the total time
    #   (e.g. since t_0 of Thornton should scale with our final end time)
    t_settled           = np.max(times) / 3e3
    valid_lums          = Luminosity[times > t_settled]
    smoothing_width     = 2 # how many checkpoints wide
    smoothing_kernel    = Gaussian1DKernel(smoothing_width)
    smoothed_lums       = convolve(valid_lums, smoothing_kernel)
    max_lum             = valid_lums[np.argmax(smoothed_lums)]
    last_run.t_0        = times[Luminosity == max_lum][0]
    last_run.t_f        = 13 * last_run.t_0 # to match with t_f given by Thornton

    parse_results = ParseResults(checkpoint_filenames)
    return parse_results

# parse_run(data_dir="../src", id="")