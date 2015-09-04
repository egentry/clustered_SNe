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
from SNe.analysis.parse import Overview



def make_dirname_from_properties(background_density, metallicity, 
                                 background_temperature, with_cooling,
                                 base_dirname=os.pardir):
    """Create a dirname where the data should live


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

def add_id_to_batch_outputs(dirname=os.pardir, batch_name=""):
    """Prepends an id to the batch output (and error) filenames, if an id exists


    Parameters
    ---------
    dirname : Optional[str]
        The name of the directory containing the batch files

    batch_name : Optional[str]
        The batch script would be named "<batch_name>.batch",
        so the output and error streams would have been captured as
        "<batch_name>*.batch.(o|e)*"

        (If batch name is an empty string, all batch files will be processed)


    Side effects
    ------------
    Renames your batch output (and error) files
    If a renamed file already exists, this *might* overwrite it 
        (See os.rename documentation for Unix/Windows behavior)


    Expects
    --------

    """

    batch_outputs = glob.glob(os.path.join(dirname, batch_name + "*.batch.o*"))
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
        dirname = make_dirname_from_properties(overview.background_density,
                                               overview.metallicity,
                                               overview.background_temperature,
                                               overview.with_cooling,
                                               base_dirname = "")
        dirname = os.path.join(target_dirname, dirname)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        files = glob.glob(os.path.join(source_dirname, overview.id + "*"))
        for file in files:
            shutil.move(file, dirname)

add_id_to_batch_outputs(dirname="../src")
# move_files(source_dirname="../tmp", target_dirname="..")




