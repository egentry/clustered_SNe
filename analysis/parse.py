import os
import glob
import warnings

import numpy as np
import pandas as pd

from astropy.convolution import convolve, Gaussian1DKernel

from scipy.signal import argrelextrema

## Boilerplate path hack to give access to full clustered_SNe package
import sys, os
if __package__ is None:
    if os.pardir not in sys.path[0]:
        file_dir = os.path.dirname(__file__)
        sys.path.insert(0, os.path.join(file_dir, 
                                        os.pardir, 
                                        os.pardir))

from clustered_SNe.analysis.constants import m_proton, pc, M_solar, \
                                             gamma, E_0, \
                                             metallicity_solar, yr

from clustered_SNe.analysis.hydro_helpers import calculate_mass, \
                                          calculate_mean_molecular_weight, \
                                          calculate_kinetic_energy, \
                                          calculate_internal_energy, \
                                          calculate_momentum, \
                                          calculate_c_ad, \
                                          calculate_entropy, \
                                          calculate_temperature, \
                                          calculate_w_cell, \
                                          calculate_crossing_time


cols = ["Radius", "dR", "dV", "Density", 
        "Pressure", "Velocity", "Z", 
        "Temperature", "Energy", "Entropy", 
        "Mass", "M_int", "C_ad", "Crossing_time"]
cols_in   = cols[:-7]


def string_to_bool(string):
    string = string.lower()
    if string in ["true", "1"]:
        return True
    if string in ["false", "0"]:
        return False

    return bool(int(string))

#####################

class RunSummary(dict):
    """Extended dict to hold useful reduced variables.

        At some point, it might be worth making this cleaner.
    """
    def __init__(self, data_dir = None, id = None):
        super(RunSummary, self).__init__()
        if data_dir is None:
            return # create an empty place holder

        checkpoint_filenames = glob.glob(os.path.join(data_dir,
                                                      id + "*checkpoint_*.dat"))
        checkpoint_filenames = sorted(checkpoint_filenames)
        num_checkpoints = len(checkpoint_filenames)

        if num_checkpoints == 0:
            raise FileNotFoundError("No checkpoints found")

        # ensure that the id is actually the FULL id
        id = get_full_id_from_partial_id(data_dir, id)

        times = np.empty(num_checkpoints)
        for k, checkpoint_filename in enumerate(checkpoint_filenames):
            f = open(checkpoint_filename, 'r')
            line = f.readline()
            times[k] = float(line.split()[3])
            f.close()

        overview = Overview(os.path.join(data_dir, id + "_overview.dat"))

        mu = calculate_mean_molecular_weight(overview.metallicity)

        E_int    = np.empty(num_checkpoints)
        E_kin    = np.empty(num_checkpoints)
        momentum = np.empty(num_checkpoints)
        E_tot    = np.empty(num_checkpoints)
        
        E_R_int  = np.empty(num_checkpoints) # Energy of the remnant
        E_R_kin  = np.empty(num_checkpoints) # Energy of the remnant
        E_R_tot  = np.empty(num_checkpoints) # Energy of the remnant
        R_shock  = np.empty(num_checkpoints) # size of the shock/remnant
        M_R      = np.empty(num_checkpoints) # mass of the shock/remnant
        
        M_tot    = np.empty(num_checkpoints)
        Z_tot    = np.empty(num_checkpoints)
        zones    = np.empty(num_checkpoints)
        
        Luminosity = np.empty(num_checkpoints)

        zones_for_each_checkpoint = [get_num_zones_in_checkpoint(checkpoint_filename)
                                     for checkpoint_filename in checkpoint_filenames]
        dataframe_columns = ["Radius",
                             "dR",
                             "dV",
                             "Density",
                             "Pressure",
                             "Velocity",
                             "Z",
                             "Mass",
                             "M_int",
                             "Temperature",
                             "Energy",
                             "Entropy",
                             "C_ad",
                             "Crossing_time",
                            ]
        tuples = [(k,i) 
                  for k,zones in enumerate(zones_for_each_checkpoint) 
                  for i in range(zones)]
        index = pd.MultiIndex.from_tuples(tuples, names=["k", "i"])
        df = pd.DataFrame(index=index, 
                          columns=dataframe_columns,
                          dtype=float)
        

        #### PARSE DATAFILES INTO DATAFRAME
        for k, checkpoint_filename in enumerate(checkpoint_filenames):
            
            df_tmp = df.loc[k]

            df_tmp.loc[:,cols_in] = pd.read_csv(checkpoint_filename,
                        delim_whitespace=True,
                        engine="c",
                        dtype=np.float64,
                        skiprows=2, 
                        names=cols_in,
                        header=None
                       ).iloc[1:-1].values

            df_tmp.Mass        = calculate_mass(df_tmp.Density.values,
                                                   df_tmp.dV.values)
            df_tmp.M_int       = df_tmp.Mass.cumsum()
            df_tmp.Temperature = calculate_temperature(df_tmp.Pressure.values,
                                                          df_tmp.Density.values, mu)
            df_tmp.Energy      = calculate_internal_energy(df_tmp.Mass.values,
                                                              df_tmp.Pressure.values,
                                                              df_tmp.Density.values) \
                                                               / df_tmp.Mass.values
            df_tmp.Entropy     = calculate_entropy(df_tmp.Temperature.values, 
                                                      df_tmp.Density.values, mu)
            df_tmp.C_ad        = calculate_c_ad(df_tmp.Pressure.values,
                                                   df_tmp.Density.values)
            w_cell = calculate_w_cell(df_tmp.Velocity.values)
            df_tmp.Crossing_time = calculate_crossing_time(df_tmp.C_ad.values,
                                                              df_tmp.Velocity.values,
                                                              w_cell,
                                                              df_tmp.dR.values )
            
            E_kin[k]    = calculate_kinetic_energy(df_tmp.Mass.values,
                                                   df_tmp.Velocity.values).sum()
            E_int[k]    = calculate_internal_energy(df_tmp.Mass.values,
                                                    df_tmp.Pressure.values,
                                                    df_tmp.Density.values).sum()
            momentum[k] = calculate_momentum(df_tmp.Mass.values,
                                             df_tmp.Velocity.values).sum()
            E_tot[k]    = E_kin[k] + E_int[k]
            M_tot[k]    = df_tmp.M_int.iloc[-1]
            Z_tot[k]    = np.sum(df_tmp.Mass * df_tmp.Z)
            zones[k]    = df_tmp.shape[0]
            
            over_dense = np.where(df_tmp.Density > overview.background_density * 1.0001)[0]
            if over_dense.size > 0:
                shock_index = over_dense.max()
            else:
                shock_index = 0

            R_shock[k] = df_tmp.iloc[shock_index].Radius
            M_R[k]     = df_tmp.iloc[shock_index].M_int

            E_R_kin[k] = calculate_kinetic_energy(df_tmp.iloc[0:shock_index+1].Mass.values,
                                                  df_tmp.iloc[0:shock_index+1].Velocity.values).sum()
            E_R_int[k] = calculate_internal_energy(df_tmp.iloc[0:shock_index+1].Mass.values, 
                                                   df_tmp.iloc[0:shock_index+1].Pressure.values,
                                                   df_tmp.iloc[0:shock_index+1].Density.values).sum()
            E_R_tot[k] = E_R_int[k] + E_R_kin[k]
            
            if k == 0:
                Luminosity[k] = 0
            else:
                Luminosity[k] = -(E_R_tot[k] - E_R_tot[k-1]) / (times[k] - times[k-1])

        df["zones"] = df.index.get_level_values(1)
        df.dR /= pc
        df.Radius /= pc
        df.Mass   /= M_solar
        df.M_int  /= M_solar
        
        self.id         = id
        self.data_dir   = data_dir
        self.df         = df
        self.times      = times
        self.zones      = zones
        self.E_tot      = E_tot
        self.Z_tot      = Z_tot
        self.E_int      = E_int
        self.E_kin      = E_kin
        self.E_R_int    = E_R_int
        self.E_R_kin    = E_R_kin
        self.E_R_tot    = E_R_tot
        self.R_shock    = R_shock
        self.M_R        = M_R
        self.momentum   = momentum
        self.Luminosity = Luminosity

        self.overview   = overview
        self.filenames  = checkpoint_filenames
        
        # filter for when initial transients have settled
        # assume that the settling time scales with the total time
        #   (e.g. since t_0 of Thornton should scale with our final end time)
        t_settled         = np.max(times) / 3e3
        valid_lums        = Luminosity[times > t_settled]
        smoothing_width   = 2 # how many checkpoints wide
        smoothing_kernel  = Gaussian1DKernel(smoothing_width)
        smoothed_lums     = convolve(valid_lums, smoothing_kernel)
        max_lum           = valid_lums[np.argmax(smoothed_lums)]
        self.t_0          = times[np.isclose(Luminosity, max_lum, atol=0)][0]
        self.t_f          = 13 * self.t_0 # to match with t_f given by Thornton

    def __repr__(self):
        """Overwrite dict.__repr__ to give the default object repr"""
        return '<%s.%s object at %s>' % (
            self.__class__.__module__,
            self.__class__.__name__,
            hex(id(self))
        )

    def __getattr__(self, name):
        return self[name]

    def __setattr__(self, name, value):
        self[name] = value

    def replace_with(self, copy_this):
        if not isinstance(copy_this, RunSummary):
            raise TypeError("Object passed to `replace_with' must be RunSummary instance")

        self.clear()
        keys = copy_this.keys()
        for key in keys:
            self[key] = copy_this[key]

    def is_last_checkpoint_x99(self):
        if self.overview.num_SNe == 0:
            return True

        last_checkpoint = self.filenames[-1]

        last_checkpoint_num = checkpoint_num_from_filename(last_checkpoint)

        return (last_checkpoint_num % 100 == 99)

    def is_converged(self):
        if self.overview.num_SNe == 0:
            return True

        momentum_max = self.momentum.max()
        tolerance = .05
        return (((1+tolerance) * self.momentum[-1]) < momentum_max)

    def is_time_resolved(self):
        if self.overview.SNe_times.size == 0:
            return True
        
        momentum_max_index = self.momentum.argmax()
        if momentum_max_index == 1:
            return False
        else:
            return True

    def is_energy_reasonable(self):
        first_unreasonable_energy_filename = self.first_unreasonable_energy()
        if first_unreasonable_energy_filename is None:
            return True
        else:
            print("unreasonable jump before checkpoint: " \
                +  first_unreasonable_energy_filename)
            return False

    def num_momentum_extrema_after_last_SNe(self, operator):
        """Operator: e.g. np.less for minima; np.greater for maxima"""
        indices = argrelextrema(self.momentum, operator)[0]
        checkpoint_after_last_SNe = np.argmax(self.times > self.overview.SNe_times.max())

        indices = indices[ indices > checkpoint_after_last_SNe]
        return len(indices)

    def first_unreasonable_energy(self):
        if self.overview.SNe_times.size == 0:
            return None

        flagged_indices = np.argwhere( (self.E_R_tot[1:]/self.E_R_tot[:-1]) > 2 )
        for i in flagged_indices:
            dE = self.E_R_tot[i+1] - self.E_R_tot[i]
            dN_SNe = sum((self.overview.SNe_times > self.times[i])
                         & (self.overview.SNe_times < self.times[i+1]))
            dE_SNe = dN_SNe * E_0
            if dE > 2*dE_SNe:
                return self.filenames[i+1]
        return None


    def num_momenta_maxima(self):
        return len(argrelextrema(self.momentum, np.greater)[0])

    def num_momenta_minima(self):
        return len(argrelextrema(self.momentum, np.less)[0])

    def momentum_peaks_too_early(self):
        if len(self.times) == 1:
            return False

        if len(self.overview.SNe_times) <= 1:
            return False

        last_SN_time = self.overview.SNe_times.max()
        last_checkpoint_time = self.times.max()
        momentum_peak_time = self.times[np.argmax(self.momentum)]

        past_last_SN        = (last_SN_time < last_checkpoint_time)
        peak_before_last_SN = (last_SN_time > momentum_peak_time)

        return (past_last_SN & peak_before_last_SN)

#####################



class Inputs(object):
    """Input Parameters"""
    def __init__(self, filename):
        super(Inputs, self).__init__()

        # read inputs properly
        self.read_inputs_from_file(filename)

    def read_inputs_from_file(self, filename):

        f = open(filename, "r")
        for line in f:
            if "T_Start" in line:
                self.T_Start = float(line.split()[-1])
            elif "T_End" in line:
                self.T_End = float(line.split()[-1])
            elif "Num_Reports" in line:
                self.Num_Reports = int(line.split()[-1])
            elif "Num_Checkpoints" in line:
                self.Num_Checkpoints = int(line.split()[-1])
            elif "Use_Logtime" in line:
                self.Use_Logtime = string_to_bool(line.split()[-1])

            elif "Num_R" in line:
                self.Num_R = int(line.split()[-1])
            elif "R_Min" in line:
                self.R_Min = float(line.split()[-1])
            elif "R_Max" in line:
                self.R_Max = float(line.split()[-1])
            elif "Log_Zoning" in line:
                self.Log_Zoning = int(line.split()[-1])
            elif "Log_Radius" in line:
                self.Log_Radius = float(line.split()[-1])

            elif "CFL" in line:
                self.CFL = float(line.split()[-1])
            elif "PLM" in line:
                self.PLM = string_to_bool(line.split()[-1])
            elif "RK2" in line:
                self.RK2 = string_to_bool(line.split()[-1])
            elif "H_0" in line:
                self.H_0 = float(line.split()[-1])
            elif "H_1" in line:
                self.H_1 = float(line.split()[-1])
            elif "Riemann_Solver" in line:
                self.Riemann_Solver = int(line.split()[-1])
            elif "Density_Floor" in line:
                self.Density_Floor = float(line.split()[-1])
            elif "Pressure_Floor" in line:
                self.Pressure_Floor = float(line.split()[-1])

            elif "With_Cooling" in line:
                self.With_Cooling = string_to_bool(line.split()[-1])
            elif "Cooling_Type" in line:
                self.Cooling_Type = line.split()[-1]
            elif "Cooling_Redshift" in line:
                self.Cooling_Redshift = float(line.split()[-1])

            elif "Adiabatic_Index" in line:
                self.Adiabatic_Index = float(line.split()[-1])

            elif "ICs" in line:
                self.ICs = line.split()[-1]

            elif "mass_loss" in line:
                self.mass_loss = line.split()[-1]

        f.close()

        if "_" in filename:
            self.id = os.path.basename(filename).split("_")[0]  
        else:
            self.id = ""

        self.dirname = os.path.dirname(filename)
        # Add trailing slash (if dirname isn't empty)
        self.dirname = os.path.join(self.dirname, "")

    def create_new_input_file(self, filename):
        f = open(filename, "w")
        f.write("T_Start: " + "{0:e}".format(self.T_Start) + "\n")
        f.write("T_End: " + "{0:e}".format(self.T_End) + "\n")
        f.write("Num_Reports: " + "{0}".format(self.Num_Reports) + "\n")
        f.write("Num_Checkpoints: " + "{0}".format(self.Num_Checkpoints) + "\n")
        f.write("Use_Logtime: " + "{0}".format(self.Use_Logtime) + "\n")

        f.write("\n")

        f.write("Num_R: " + "{0}".format(self.Num_R) + "\n")
        f.write("R_Min: " + "{0:e}".format(self.R_Min) + "\n")
        f.write("R_Max: " + "{0:e}".format(self.R_Max) + "\n")
        f.write("Log_Zoning: " + "{0}".format(self.Log_Zoning) + "\n")
        f.write("Log_Radius: " + "{0}".format(self.Log_Radius) + "\n")

        f.write("\n")


        f.write("CFL: " + "{0}".format(self.CFL) + "\n")
        f.write("PLM: " + "{0}".format(self.PLM) + "\n")
        f.write("RK2: " + "{0}".format(self.RK2) + "\n")
        f.write("H_0: " + "{0}".format(self.H_0) + "\n")
        f.write("H_1: " + "{0}".format(self.H_1) + "\n")
        f.write("Riemann_Solver: " + "{0}".format(self.Riemann_Solver) + "\n")
        f.write("Density_Floor: " + "{0}".format(self.Density_Floor) + "\n")
        f.write("Pressure_Floor: " + "{0}".format(self.Pressure_Floor) + "\n")

        f.write("\n")

        f.write("With_Cooling: " + "{0}".format(self.With_Cooling) + "\n")
        f.write("Cooling_Type: " + "{0}".format(self.Cooling_Type) + "\n")
        f.write("Cooling_Redshift: " + "{0}".format(self.Cooling_Redshift) + "\n")

        f.write("\n")

        f.write("Adiabatic_Index: " + "{0}".format(self.Adiabatic_Index) + "\n")

        f.write("\n")

        f.write("ICs: " + "{0}".format(self.ICs) + "\n")

        f.write("\n")

        f.write("mass_loss: " + "{0}".format(self.mass_loss) + "\n")


#####################


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
        """Create an Overview object using an "*overview.dat" filename

        Parameters
        ----------
        filename : str
            - should be a valid "overview.dat" style filename (see Output/ascii.cxx)

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
            Cooling Type : str
            Number of SNe : int
            Cluster Mass : float
                [g]
            seed : int
            Mass Loss : str

        """
        super(Overview, self).__init__()
        
        if "_" in filename:
            self.id = os.path.basename(filename).split("_")[0]  
        else:
            self.id = ""
        self.dirname = os.path.dirname(filename)
        # Add trailing slash (if dirname isn't empty)
        self.dirname = os.path.join(self.dirname, "")
        
        inputs_filename = os.path.join(os.path.dirname( filename),
                                       os.path.basename(filename).replace("overview", "inputs"))
        self.inputs = Inputs(inputs_filename)

         # default, since earlier runs won't have this saved
        self.cluster_mass = 0
        self.cooling_type = "equilibrium"
        self.mass_loss = "none"

        f = open(filename, "r")
        for line in f:
            if "Metallicity" in line:
                self.metallicity = float(line.split()[1])
            elif "Background Density" in line:
                self.background_density = float(line.split()[2])
            elif "Background Temperature" in line:
                self.background_temperature = float(line.split()[2])
            elif "With cooling" in line:
                self.with_cooling = string_to_bool(line.split()[2])
            elif "Cooling Type" in line:
                self.cooling_type = line.split()[-1]
            elif "Number of SNe" in line:
                self.num_SNe = int(line.split()[-1])
            elif "Cluster Mass" in line:
                self.cluster_mass = float(line.split()[-1]) * M_solar
            elif "seed" in line:
                self.seed = int(line.split()[-1])
            elif "mass loss" in line: 
                self.mass_loss = line.split()[-1]
        f.close()

        SNe_filename = os.path.join(os.path.dirname( filename),
                                    os.path.basename(filename).replace("overview", "SNe"))
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='loadtxt: Empty input file: "*"')
            SNe = np.loadtxt(SNe_filename, ndmin=2, usecols=range(5))
        if (SNe.shape[0] != self.num_SNe):
            raise ValueError("Number of SNe in datafile " +
                             "doesn't match number listed in overview file" +
                             " for file: " + filename)
        if self.num_SNe > 0:
            self.SNe_times          = SNe[:,0]
            self.SNe_initial_mass   = SNe[:,1]
            self.SNe_ejecta_mass    = SNe[:,2]
            self.SNe_ejecta_mass_Z  = SNe[:,3]
            self.SNe_wind_mass      = SNe[:,4]

            sorted_indices = np.argsort(self.SNe_times)
            self.SNe_times          = self.SNe_times[         sorted_indices]
            self.SNe_initial_mass   = self.SNe_initial_mass[  sorted_indices]
            self.SNe_ejecta_mass    = self.SNe_ejecta_mass[   sorted_indices]
            self.SNe_ejecta_mass_Z  = self.SNe_ejecta_mass_Z[ sorted_indices]
            self.SNe_wind_mass      = self.SNe_wind_mass[     sorted_indices]
        else:
            self.SNe_times          = np.array([])
            self.SNe_initial_mass   = np.array([])
            self.SNe_ejecta_mass    = np.array([])
            self.SNe_ejecta_mass_Z  = np.array([])
            self.SNe_wind_mass      = np.array([])

        return

    def __str__(self):
        string  = "id \t\t\t = {0}".format(self.id) + "\n"
        string += "metallicity \t\t = {0} ".format(self.metallicity) + "\n"
        string += "background density \t = {0} [g cm^-3]".format(self.background_density) + "\n"
        string += "background temperature \t = {0:e} [K]".format(self.background_temperature)
        return string


#####################



def extract_masses_momenta_raw(data_dir, density, metallicity,
                           H_1=.1,
                           extract_at_last_SN = False):
    """For a faster version, get these from the database,
    using `extract_masses_momenta` from database_helpers"""
    
    if not os.path.exists(data_dir):
        raise FileNotFoundError("No directory found named: "+ data_dir)
        
    overview_filenames = glob.glob(os.path.join(data_dir,
                                                "*overview.dat"))
    ids = np.array([])
    masses = np.array([])
    momenta = np.array([])
    
    
    for overview_filename in overview_filenames:
        overview = Overview(overview_filename)
        if extract_at_last_SN and (overview.num_SNe < 2):
            continue
        
        
        # apply filters
        
        if overview.inputs.H_1 != H_1:
            continue
        
        if not np.isclose(overview.background_density, density, atol=0):
            continue
            
        if not np.isclose(overview.metallicity, metallicity, atol=0):
            continue
        
        mass = overview.cluster_mass
        id = os.path.basename(overview_filename).split("_")[0]
    
        run_summary = RunSummary(data_dir=data_dir, id=id)
        if extract_at_last_SN:
            last_SNe_time = run_summary.overview.SNe_times.max()
            times = run_summary.times
            
            checkpoint_after_SNe = np.argmin(np.abs(times - times[times>last_SNe_time].min()))
            momentum = run_summary.momentum[checkpoint_after_SNe]
        else:
            momentum = run_summary.momentum.max()
        
        masses  = np.append(masses, mass)
        momenta = np.append(momenta, momentum)
        ids     = np.append(ids, id)

    if len(ids) == 0:
        print("No matching files found!")
        return
    
    sorted_indices = np.argsort(masses)
    masses  = masses[  sorted_indices]
    momenta = momenta[ sorted_indices]
    ids     = ids[     sorted_indices]
    return masses, momenta, ids

def get_full_id_from_partial_id(data_dir, partial_id):
    filenames = glob.glob(os.path.join(data_dir, partial_id + "*_checkpoint_*.dat"))
    if len(filenames) == 0:
        raise FileNotFoundError("No checkpoints found with partial id: " + partial_id)

    full_id = os.path.basename(filenames[0]).split("_checkpoint")[0]

    for filename in filenames:
        full_id_tmp = os.path.basename(filename).split("_checkpoint")[0]
        if full_id_tmp != full_id:
            raise ValueError("partial id not unique in directory")
    return full_id

def checkpoint_num_from_filename(checkpoint_filename):
    if "checkpoint" not in checkpoint_filename:
        raise ValueError("checkpoint filename '" + checkpoint_filename + "' not valid checkpoint file")
    return int(checkpoint_filename.split("_")[-1].strip(".dat"))


def get_num_zones_in_checkpoint(checkpoint_filename):
    num_header_lines = 2
    num_guard_cells = 2
    with open(checkpoint_filename) as f:
        for i, l in enumerate(f):
            pass
    return i - num_header_lines - num_guard_cells + 1


def parse_into_scientific_notation(number, prefix_format="{:.2e}"):
    """Doesn't include the bounding $ for the math mode"""
    scientific_notation_parts = prefix_format.format(number).split("e")
    scientific_notation_parts[1] = int(scientific_notation_parts[1])

    return r"{0} \cdot 10^{{{1}}}".format(scientific_notation_parts[0],
                                         scientific_notation_parts[1])

