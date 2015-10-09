
from __future__ import print_function, division

import os, sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact, interactive, interact_manual
from ipywidgets import FloatSlider, IntSlider, Dropdown
from ipywidgets import Checkbox, RadioButtons, fixed 

from IPython.display import display, Math, clear_output


from astropy.convolution import convolve, Gaussian1DKernel

import seaborn as sns


## import from local files
## Boilerplate path hack to give access to full SNe package
import sys
if __package__ is None:
    if os.pardir not in sys.path[0]:
        file_dir = os.getcwd()
        sys.path.insert(0, os.path.join(file_dir, 
                                        os.pardir, 
                                        os.pardir))

from SNe.analysis.constants import m_proton, pc, yr, M_solar, \
                                   gamma, E_0, metallicity_solar
    
from SNe.analysis.sedov.dimensionalize_sedov import dimensionalized_sedov
from SNe.analysis.sedov.closed_form_sedov import SedovSolution
from SNe.analysis.parameter_study_file_structure import make_dirname_from_properties
from SNe.analysis.parse import RunSummary, Overview, parse_run, cols


general_string_format = ".2e"



cols_plot = cols 
cols_plot_linear = [cols[i] for i in [6, 9]]

label_dict = {"Radius"       : r"$R$ [pc]", 
              "Velocity"     : r"$U$ [cm s$^{-1}$]",
              "Density"      : r"$\rho$  [g cm$^{-3}$]",
              "Temperature"  : r"$T$ [K]",
              "Mass"         : r"$M$ [$M_\odot$]",
              "M_int"        : r"$M_\mathrm{int}$ [$M_\odot$]",
              "C_ad"         : r"$C_{ad}$ [cm s$^{-1}$]",
              "Crossing_time": r"$dt_{\mathrm{crossing}}$",
              "Energy"       : r"$E_{int}$ [erg g$^{-1}$]",
              "Pressure"     : r"$P$ [dyne cm$^{-2}$]",
              "Entropy"      : r"$S$ [$k_B$ / particle]",
              "M_int"        : r"$M_{int}$ [$M_\odot$]",
              "dR"           : r"$\Delta R$ [cm]",
              "dV"           : r"$\Delta V$ [cm$^3$]",
              "Z"            : r"$Z$ (metallicity)",
              "zones"        : r"zones"
             }



def plotter(last_run,
            checkpoint_filenames, 
            metallicity, background_density, background_temperature,
            sedov_solution,
            x_axis_variable  = "Radius",
            y_axis_variable  = "Density",
            with_Sedov       = True,
            highlight_timestep_limiting_cell = False,
            outer_limit_log  = 0, 
            checkpoint_index = 0):
    df_tmp = last_run["df"].loc[checkpoint_index]

    checkpoint_filename = checkpoint_filenames[checkpoint_index]
    time = last_run["times"][checkpoint_index]
    
    E_kin = sedov_solution.E_kin
    E_int = sedov_solution.E_int
    print("E_kin: ", format(E_kin, general_string_format))
    print("E_int: ", format(E_int, general_string_format))
    momentum = sedov_solution.get_momentum(time=time - last_run["times"][0])

    print("checkpoint: ",
          checkpoint_filename)
    print("time:                      ",
        format(time / yr, general_string_format), "[yr]",
        "\t", format(time, general_string_format), "[s]")
    print("time elapsed:              ",
          format((time - last_run["times"][0]) / yr, general_string_format), "[yr]",
          "\t", format((time - last_run["times"][0]), general_string_format), "[s]")

    if last_run["overview"].SNe_times.size == 1:
        print("energy conserved to:       ", 
              format( (   last_run["E_tot"][checkpoint_index]
                        - last_run["E_tot"][0])
                      / last_run["E_tot"][0], general_string_format) )
        print("E_kin    accurate to:      ", 
              format( (last_run["E_kin"][checkpoint_index]
                              - E_kin)
                      / E_kin, general_string_format) )
        print("momentum accurate to:      ", 
              format( (last_run["momentum"][checkpoint_index]
                             - momentum)
                      / momentum, general_string_format) )
        # # DON'T USE E_INT AS A METRIC, since you accrete in thermal energy
    #     print("E_int accurate to:         ", 
    #           (last_run["E_int"][checkpoint_index] - E_int) / E_int)
        print("Peak luminosity at checkpoint",
              np.argmax(last_run["times"] == last_run["t_0"]) )
        print("Peak luminosity at t_0 =   ",
              format(last_run["t_0"] / yr, general_string_format), "[yr]")
        print("t_f = 13 * t_0 =           ",
              format(last_run["t_f"] / yr, general_string_format), "[yr]")
        print("R_shock =                  ",
              format(last_run["R_shock"][checkpoint_index] / pc, "3.2f"), "[pc]")
        print("E_R_tot =                  ",
              format(last_run["E_R_tot"][checkpoint_index], general_string_format),
              "[ergs]")
    print("background_density:        ", 
        format(last_run["overview"].background_density, general_string_format))
    print("Cluster mass:              ", 
        format(last_run["overview"].cluster_mass / M_solar, general_string_format),
        "M_sol")
    print("Number of SNe so far:      ",
          np.sum(last_run["overview"].SNe_times <= time))
    print("mass loss prescription:    ", last_run["overview"].mass_loss)
    
    if x_axis_variable is "Radius":
        plt.xlim((0,10**outer_limit_log))
    if x_axis_variable is "M_int":
        plt.xlim((0,10**2))

    marker = "."
    
    plt.plot(df_tmp[x_axis_variable], df_tmp[y_axis_variable], 
             marker=marker,
             label="numeric",
             drawstyle="steps")
    if highlight_timestep_limiting_cell is True:
        timestep_limiting_index = df_tmp.Crossing_time.argmin()
        plt.plot(df_tmp[x_axis_variable].loc[timestep_limiting_index],
                 df_tmp[y_axis_variable].loc[timestep_limiting_index], 
                 marker=marker,
                 linestyle="",
                 label="Timestep limiting cell",
                 color='r')

    if y_axis_variable in cols_plot_linear:
        plt.yscale("linear")
    else:
        plt.yscale("log")
        # Set and fix limits
        y_min, y_max = plt.ylim()
        plt.ylim( (y_min / 5, y_max * 5) )
    if y_axis_variable is "Velocity":
        plt.ylim(ymin=1)
    
    if with_Sedov:
        plot_sedov(last_run, time, x_axis_variable, y_axis_variable, 
                   metallicity, background_density, background_temperature)
    
    plt.xlabel(label_dict[x_axis_variable])
    plt.ylabel(label_dict[y_axis_variable])

    plt.legend(loc="best")
    

def plot_sedov(last_run, time, x_axis_variable, y_axis_variable, 
               metallicity, background_density, background_temperature):

    sedov_x_axes = ["Radius", "M_int"]
    if x_axis_variable not in sedov_x_axes:
        return
    
    sedov_cols = ["Radius", "Velocity", "Density", "Temperature",
                  "C_ad", "Energy", "Pressure", "Entropy", "Mass"]
    sedov_cols_plot = sedov_cols[:-1]
    if y_axis_variable not in sedov_cols_plot:
        return
    
    sedov_time = time - last_run["overview"].SNe_times[0]
    if sedov_time <= 0:
        return
    
    SNe_so_far = np.sum(last_run["overview"].SNe_times <= time) 
    if SNe_so_far != 1:
        return
    
    sedov = dimensionalized_sedov(time - last_run["times"][0],
                                  metallicity=metallicity, 
                                  background_density=background_density,
                                  background_temperature=background_temperature)
    sedov = np.array(sedov).transpose()
    df_sedov = pd.DataFrame.from_records(sedov, 
                                         columns=sedov_cols)
    df_sedov.Radius  /= pc
    df_sedov.Mass    /= M_solar
    df_sedov["M_int"] = df_sedov.Mass.cumsum()

    plt.plot(df_sedov[x_axis_variable], df_sedov[y_axis_variable], 
             label="analytic (no cooling)")
    
def single_run(data_dir="", id=""):

    if not os.path.exists(data_dir):
        raise FileNotFoundError("No directory found named: "+ data_dir)
        
    last_run, parse_results = parse_run(data_dir, id)
    sedov_solution = SedovSolution(E_0,
                                   last_run["overview"].background_density, 
                                   last_run["overview"].metallicity)
    
    #### PASS TO PLOTTER ####
    num_checkpoints = len(parse_results.checkpoint_filenames)
    
    log_R_max = round(np.log10(last_run["df"]["Radius"].max()), 2)
    log_R_min = max(log_R_max-4, 
                    round(np.log10(last_run["df"]["Radius"].min()), 2)+1)
                
    if type(single_run.previous_widget) is widgets.Box:
        single_run.previous_widget.close()

    w = interactive(plotter,
        last_run               = fixed(last_run),
        checkpoint_filenames   = fixed(parse_results.checkpoint_filenames),
        metallicity            = fixed(last_run["overview"].metallicity),
        background_density     = fixed(last_run["overview"].background_density),
        background_temperature = fixed(last_run["overview"].background_temperature),
        sedov_solution         = fixed(sedov_solution),
        outer_limit_log        = FloatSlider(min=log_R_min, 
                                             max=log_R_max, 
                                             step=0.1, 
                                             value=log_R_max),
        checkpoint_index       = IntSlider(min=0, 
                                           max=num_checkpoints-1, 
                                           step=0, 
                                           value=num_checkpoints-1),
        y_axis_variable        = Dropdown(options=cols_plot,
                                          value="Density"),
        x_axis_variable        = RadioButtons(options=["Radius",
                                                       "M_int",
                                                       "zones"]))
    single_run.previous_widget = w
    display(w)
    return last_run
single_run.previous_widget = widgets.Box()

def conduction_comparisons(mass, H_0, data_dir,
							num_SNe=1, 
							plot_these_H_1=[0, .1, .3, 1, 10]):
        
    if not os.path.exists(data_dir):
        raise FileNotFoundError("No directory found named: "+ data_dir)
        
    overview_filenames = glob.glob(os.path.join(data_dir,
                                                "*overview.dat"))
    
    if H_0 is True:
        label = r"$H_0"
    else:
        label = r"$H_1"

    ids = np.array([])
    labels = np.array([])
    conduction_values = np.array([])
    for overview_filename in overview_filenames:
        overview = Overview(overview_filename)
        if overview.cluster_mass != mass:
            continue
        if mass == 100 * M_solar:
            if overview.num_SNe != num_SNe:
                continue
            
        ids = np.append(ids, 
                        os.path.basename(overview_filename).split("_")[0])
        if H_0 is True:
            labels = np.append(labels,
                               label + r" = {0}$".format(int(round(overview.inputs.H_0*10)) / 10))
            conduction_values = np.append(conduction_values,
                                          overview.inputs.H_0)
        else:
            if overview.inputs.H_1 not in plot_these_H_1:
                continue
            labels = np.append(labels,
                               label + r" = {0}$".format(int(round(overview.inputs.H_1*10)) / 10))
            conduction_values = np.append(conduction_values,
                                          overview.inputs.H_1)

    if len(ids) == 0:
        print("No matching files found!")
        return
            
    sorted_indices = np.argsort(conduction_values)
    conduction_values = conduction_values[sorted_indices]
    ids = ids[sorted_indices]
    labels = labels[sorted_indices]

    last_checkpoints = np.array([], dtype=int)
    last_common_checkpoint = np.inf
    for i, id in enumerate(ids):
        num_checkpoint = len(glob.glob(os.path.join(data_dir,
                                                    id+"*checkpoint*")))
        if num_checkpoint < last_common_checkpoint:
            last_common_checkpoint = num_checkpoint
        last_checkpoints = np.append(last_checkpoints,
                                     num_checkpoint)

    for i, id in enumerate(ids):
        print("====")
        display(Math(labels[i]))
        print("Number of checkpoints available: ",
              last_checkpoints[i])
        last_run, parse_results = parse_run(data_dir, id)
        sedov_solution = SedovSolution(E_0,
                                       last_run["overview"].background_density, 
                                       last_run["overview"].metallicity)

        #### PASS TO PLOTTER ####
        num_checkpoints = len(parse_results.checkpoint_filenames)
        plot_checkpoint = last_common_checkpoint-1

        log_R_max = round(np.log10(last_run["df"]["Radius"].max()), 2)
        log_R_min = max(log_R_max-4, 
                        round(np.log10(last_run["df"]["Radius"].min()), 2)+1)


        if last_run["overview"].num_SNe == 1:
            SN_or_SNe = "SN"
        else:
            SN_or_SNe = "SNe"
        plt.title("Num " + SN_or_SNe + ": {0}".format(last_run["overview"].num_SNe))

        plotter(last_run,
                parse_results.checkpoint_filenames, 
                last_run["overview"].metallicity, 
                last_run["overview"].background_density, 
                last_run["overview"].background_temperature,
                sedov_solution,
                x_axis_variable  = "Radius",
                y_axis_variable  = "Temperature",
                with_Sedov       = False,
                highlight_timestep_limiting_cell = True,
                outer_limit_log  = log_R_max, 
                checkpoint_index = plot_checkpoint)
    
    ax = plt.gca()
    handles_tmp, labels_tmp = ax.get_legend_handles_labels()
    handles = handles_tmp[::2]
    handles = np.append(handles, 
                        handles_tmp[-1])
    labels = np.append(labels,
                       "Timestep limiting cells")
    ax.legend(handles, labels, loc="best")

    plot_filename = "plots/conduction"
    plot_filename += "_{0:d}_M_sun".format(int(round(mass / M_solar)))
    plot_filename += "_{0:s}".format(label.strip("$"))
    print("plot_filename:", plot_filename)
    plt.savefig(plot_filename + ".eps")
    plt.savefig(plot_filename + ".pdf")
    plt.savefig(plot_filename + ".png")


def parameter_study_wrapper(log_n, log_Z, T=1e4, 
                            with_cooling=True):
    z_solar = metallicity_solar
    
    background_density     = 1.33 * m_proton * 10**log_n
    metallicity            = z_solar * 10**log_Z
    background_temperature = T
    
    data_dir = make_dirname_from_properties( background_density, 
                                             metallicity, 
                                             background_temperature,
                                             with_cooling)
    last_run = RunSummary() # default
    if os.path.isdir(data_dir):
        tmp = glob.glob(os.path.join(data_dir, "*checkpoint*.dat"))        
        
        if len(tmp) is not 0:
            basename = os.path.basename(tmp[0])
            id = basename.split("checkpoint")[0]
            print(data_dir)
            last_run = single_run(data_dir=data_dir, id=id)

                        
        else:
            print("No data was found")
    else:
        print("No directory was found")
        print("Directory: ", data_dir)
    return last_run




def SNe_distplot(last_run, x_axis):
    if x_axis is "time":
        x_data = last_run["overview"].SNe_times / yr
        rug=True
        hist=False
    elif x_axis is "checkpoints":
        x_data = np.array([], dtype=np.int)
        
        for SNe_time in last_run["overview"].SNe_times:
            if (SNe_time >= last_run["times"].min()) and (SNe_time < last_run["times"].max()):
                x_data = np.append(x_data, np.argmin(np.abs(last_run["times"] - SNe_time)))
        # this would be more natural as a histogram,
        # but I can't figure out how to normalize a histogram in a good way
        # Maybe it'd just be better to use subplots?
        rug=True
        hist=False
    else:
        raise ValueError("Unrecognized value for x_axis: " + x_axis )
    
    if x_data.size == 1:
        x_data = np.tile(x_data, 2) #seaborn can't do a rug plot of 1 point
    sns.distplot(x_data, color="k", norm_hist=False, 
                 hist=hist, rug=rug, kde=False, 
                 rug_kws={"linewidth":3},
                 bins=np.arange(last_run["times"].size))


def plot_zones(last_run, distplot=True):
    if type(last_run) is RunSummary:
        if last_run["zones"] is not None:
            plt.plot(last_run["zones"])
            plt.ylim(ymin=0)
            plt.ylabel("Number of Zones")
            plt.xlabel("Checkpoint")
            if distplot is True:
                SNe_distplot(last_run, "checkpoints")

def plot_shock_location(last_run, clear_previous = True, distplot=True):
    if clear_previous:
        plt.figure()
    if distplot is True:
        SNe_distplot(last_run, "time")
    plt.plot(last_run["times"] / yr, last_run["R_shock"] / pc)
    plt.xlabel(r"time [yr]")
    plt.ylabel(r"$R_{\mathrm{shock}}$ [pc]")

def plot_energy(last_run, x_axis):
    if type(last_run) is RunSummary:
        if ( (last_run["E_tot"] is not None) and
             (last_run["E_int"] is not None) and
             (last_run["E_kin"] is not None) ):
            
            plt.figure()
            
            if x_axis is "time":
                x_variable = last_run["times"] / yr
                xlabel = "Time [yr]"
                xscale = "linear"
                plt.xscale(xscale)
                xfmt = plt.gca().get_xaxis().get_major_formatter() # needs to be set AFTER plt.xscale()
                if xscale is "log":
                    mask = x_variable > 1
                elif xscale is "linear":
                    mask = np.full_like(x_variable, True, dtype=bool) 
                    xfmt.set_powerlimits((-2, 2)) # force scientific notation outside this range

            else:
                x_variable = np.arange(len(last_run["times"]))
                xlabel = "Checkpoint"
                xscale = "linear"
                mask = np.full_like(x_variable, True, dtype=bool) 

                plt.xscale(xscale)
                xfmt = plt.gca().get_xaxis().get_major_formatter() # needs to be set AFTER plt.xscale()
                
                
            
            E_err = (last_run["E_tot"] - last_run["E_tot"][0]) / last_run["E_tot"][0]
            plt.plot(x_variable[mask], E_err[mask])
            plt.xscale(xscale)
            plt.xlabel(xlabel)   
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.ylabel("Fractional Change (Energy)")
            SNe_distplot(last_run, x_axis)

            plt.figure()
            plt.plot(x_variable[mask], last_run["E_tot"][mask], label="E_tot" )
            plt.plot(x_variable[mask], last_run["E_kin"][mask], label="E_kin" )
            plt.plot(x_variable[mask], last_run["E_int"][mask], label="E_int" )
            plt.legend(loc="best")
            plt.xscale(xscale)
            plt.xlabel(xlabel) 
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.ylabel("Energy [erg]")
            SNe_distplot(last_run, x_axis)

            
            plt.figure()
            plt.plot(x_variable[mask], last_run["E_R_tot"][mask], label="E_Remnant" )
            plt.legend(loc="best")
            plt.xscale(xscale)
            plt.xlabel(xlabel)  
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.ylabel("Energy [erg]")
            SNe_distplot(last_run, x_axis)


def plot_momentum(last_run, x_axis, clear_previous=True, distplot=True):
    if last_run["overview"].cluster_mass <= 0:
        raise ValueError("Cluster mass doesn't allow valid normalization of momentum")
    if type(last_run) is RunSummary:
        if clear_previous:
            plt.figure()

        if x_axis is "time":
            x_variable = last_run["times"] / yr
            xlabel = "Time [yr]"
            xscale = "linear"
            plt.xscale(xscale)
            xfmt = plt.gca().get_xaxis().get_major_formatter() # needs to be set AFTER plt.xscale()
            if xscale is "log":
                mask = x_variable > 1
            elif xscale is "linear":
                mask = np.full_like(x_variable, True, dtype=bool) 
                xfmt.set_powerlimits((-2, 2)) # force scientific notation outside this range

        else:
            x_variable = np.arange(len(last_run["times"]))
            xlabel = "Checkpoint"
            xscale = "linear"
            mask = np.full_like(x_variable, True, dtype=bool) 

            plt.xscale(xscale)
            
            # needs to be set AFTER plt.xscale():
            xfmt = plt.gca().get_xaxis().get_major_formatter() 
        
        if distplot is True:
            SNe_distplot(last_run, x_axis)


        plt.plot(x_variable[mask], 
                 last_run["momentum"][mask] / (last_run["overview"].cluster_mass * 100*1000))
        plt.xscale(xscale)
        plt.xlabel(xlabel)   
        plt.gca().xaxis.set_major_formatter(xfmt)
        plt.ylabel(r"Momentum / M$_\mathrm{cluster}$ [km s$^{-1}$]")
        
        plt.ylim(ymin=0)



def plot_luminosity(last_run, x_axis):
    if type(last_run) is RunSummary:
        plt.figure()
        if x_axis is "time":
            x_variable = last_run["times"] / yr
            xlabel = "Time [yr]"
            xscale = "log"
            if xscale is "log":
                mask = x_variable > 1
            else:
                mask = np.full_like(x_variable, True, dtype=bool)
                xfmt.set_powerlimits((-2, 2)) # force scientific notation outside this range

        else:
            x_variable = np.arange(len(last_run["times"]))
            xlabel = "Checkpoint"
            xscale = "linear"
            mask = np.full_like(x_variable, True, dtype=bool)
        
        y_data = last_run["Luminosity"][mask]
        gauss_kernel = Gaussian1DKernel(2)
        y_data = convolve(y_data, gauss_kernel)
        
        plt.plot(x_variable[mask],  y_data, label="net cooling")
        plt.plot(x_variable[mask], -y_data, label="net heating + accretion", color="r")
        plt.xlabel(xlabel)
        plt.xscale(xscale)
        plt.ylabel(r"dE/dt [erg s$^{-1}$]")
        plt.yscale("log")
        plt.legend(loc="best")
        
        xfmt = plt.gca().get_xaxis().get_major_formatter()

        plt.gca().xaxis.set_major_formatter(xfmt)
        print("Luminosity max at checkpoint: ", np.argmax(last_run["t_0"] == last_run["times"]))
        print("Luminosity max at time:       ", format(last_run["t_0"] / yr,
                                                       general_string_format),
              "[yr]" )