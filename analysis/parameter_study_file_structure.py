import os
import glob
import shutil
import numpy as np

## Boilerplate path hack to give access to full SNe package
import sys
if __package__ is None:
    if os.pardir not in sys.path[0]:
        file_dir = os.path.dirname(__file__)
        sys.path.insert(0, os.path.join(file_dir, 
                                        os.pardir, 
                                        os.pardir))

from SNe.analysis.constants import m_proton

class Overview(object):
    """Overview of a given ./SNe run


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

        """
        super(Overview, self).__init__()
        
        self.id = os.path.basename(filename).split("_")[0]  
        self.dirname = os.path.dirname(filename)
        # Add trailing slash (if dirname isn't empty)
        self.dirname = os.path.join(self.dirname, "")
        
        self.num_SNe = 0 # default, since earlier runs won't have this saved

        # this parsing is redundant with parts from parse_run in parse.py
        # they should probably be consolidated
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
        f.close()

        SNe_times_filename = filename.replace("overview", "SNe_times") 
        if os.path.exists(SNe_times_filename):
            self.SNe_times = np.loadtxt(SNe_times_filename)
        else:
            self.SNe_times = np.array([0.])
        self.SNe_times.sort()


        return

    def __str__(self):
        string  = "id \t\t\t = {0}".format(self.id) + "\n"
        string += "metallicity \t\t = {0} ".format(self.metallicity) + "\n"
        string += "background density \t = {0} [g cm^-3]".format(self.background_density) + "\n"
        string += "background temperature \t = {0:e} [K]".format(self.background_temperature)
        return string

    def make_dirname(self):
        """Create a dirname where the data for this Overview should live


        Expects
        --------
        The following attributes should be set:
            dirname 
            metallicity
            background_density
            background_temperature
            with_cooling


        Returns
        -------
        dirname : str
            dirname is not required to exist


        Side effects
        ------------
        None

        """

        dirname = make_dirname_from_properties(self.background_density, 
                                               self.metallicity, 
                                               self.background_temperature,
                                               self.with_cooling,
                                               base_dirname = "")
        return dirname


def make_dirname_from_properties(background_density, metallicity, 
                                 background_temperature, with_cooling,
                                 base_dirname=os.pardir):
    """Create a dirname where the data for this Overview should live


    Parameters
    --------
        background_density : float
        metallicity : float
        background_temperature : float
        with_cooling : bool
        base_dirname : Optional[str]
            A properly constructed dirname, on which we build a relative dirname 
            for a specific run


    Returns
    -------
    dirname : str
        This directory is not required to exist


    Side effects
    ------------
    None

    """

    dirname = os.path.join(base_dirname, "saved_runs", "riemann", "Thornton_parameter_study")
    dirname = os.path.join(dirname, "log_n_" + "{0:+03.0f}".format(5*int(round(2*np.log10(background_density / m_proton / 1.33)))))
    dirname = os.path.join(dirname, "log_Z_" + "{0:+03.0f}".format(5*int(round(2*np.log10(metallicity / .02)))))
    dirname = os.path.join(dirname, "T_" + "{0:.0f}".format(background_temperature))
    if with_cooling is True:
        dirname = os.path.join(dirname, "with_cooling")
    else:
        dirname = os.path.join(dirname, "no_cooling")
    # Add trailing slash
    dirname = os.path.join(dirname, "")

    return dirname

def add_id_to_batch_outputs(dirname=os.pardir, batch_name="Thornton_parameter_study"):
    """Prepends an id to the batch output (and error) filenames, if an id exists


    Parameters
    ---------
    dirname : Optional[str]
        The name of the directory containing the batch files

    batch_name : Optional[str]
        The batch script would be named "<batch_name>.batch",
        so the output and error streams would have been captured as
        "<batch_name>.batch.(o|e)*"


    Side effects
    ------------
    Renames your batch output (and error) files
    If a renamed file already exists, this *might* overwrite it 
        (See os.rename documentation for Unix/Windows behavior)


    Expects
    --------

    """

    batch_outputs = glob.glob(os.path.join(dirname, batch_name + ".batch.o*"))
    for batch_output in batch_outputs:
        f = open(batch_output, "r")
        for line in f:
            if "uuid" in line:
                id = line.split()[-1]
                if id not in batch_output:
                    batch_output_parts = os.path.split(batch_output)
                    new_batch_output = os.path.join(batch_output_parts[0],
                                             id + "_" + batch_output_parts[1])
                    os.rename(batch_output, new_batch_output)
                batch_error = batch_output.replace("batch.o", "batch.e")
                if id not in batch_error:
                    batch_error_parts = os.path.split(batch_error)
                    new_batch_error = os.path.join(batch_error_parts[0],
                                             id + "_" + batch_error_parts[1])
                    os.rename(batch_error, new_batch_error)
                break
        f.close()

    return


def move_files(source_dirname=os.path.join(os.pardir, "src"),
               target_dirname=os.path.join(os.pardir)):
    """Moves files from current location into their locations on a tree of directories


    Parameters
    ---------
    source_dirname : Optional[str]
        The name the directory containing the files initially

    target_dirname: Optional[srt]
        Files will be saved to "<target_dirname>/saved_runs/riemann..."


    Side effects
    ------------
    Renames your batch output (and error) files
    If a renamed file already exists, this *might* overwrite it 
        (See os.rename documentation for Unix/Windows behavior)


    Expects
    --------
    Output and error files should start with "<batch_name>.batch.(o|e)"
        This is the default behavior if you submit a "<batch_name>.batch" file to SGE qsub

    """
    add_id_to_batch_outputs(dirname=source_dirname)
    overview_filenames = glob.glob(os.path.join(source_dirname, "*overview.dat"))
    overviews = [None] * len(overview_filenames)
    for i, overview_filename in enumerate(overview_filenames):
        overviews[i] = Overview(overview_filename)

    for overview in overviews:
        dirname = overview.make_dirname()
        dirname = os.path.join(target_dirname, dirname)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        files = glob.glob(os.path.join(source_dirname, overview.id + "*"))
        for file in files:
            shutil.move(file, dirname)

# move_files(source_dirname="../src", target_dirname="..")



