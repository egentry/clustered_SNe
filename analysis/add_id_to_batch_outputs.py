import sys, os
import glob

def add_id_to_batch_outputs(dirname="", batch_name=""):
    """Prepends an id to the batch output (and error) filenames, if an id exists

    Parameters
    ---------
    dirname : Optional[str]
        The name of the directory containing the batch files

    batch_name : Optional[str]
        The batch script would be named "<batch_name>.batch",
        so the output and error streams would have been captured as
        "<batch_name>*.batch*.(o|e)*"

        (If batch name is an empty string, all batch files will be processed)


    Side effects
    ------------
    Renames your batch output (and error) files
    If a renamed file already exists, this *might* overwrite it 
        (See os.rename documentation for Unix/Windows behavior)


    """
    batch_outputs = glob.glob(os.path.join(dirname, batch_name + "*.batch*.o*"))

    for batch_output in batch_outputs:
        f = open(batch_output, "r")
        for line in f:
            if ("uuid" in line) or ("output_prefix" in line):
                id = line.split()[-1].strip("_")
                id += "_"
                
                dir, batch_output_base = os.path.split(batch_output)

                #strip id, just in case
                batch_output_base = batch_output_base.replace(id, "")

                batch_output_new = os.path.join(dir, id + batch_output_base)

                try: 
                    os.rename(batch_output, batch_output_new)
                except FileNotFoundError:
                    pass

                # this isn't terribly robust
                batch_error_base     = batch_output_base.replace(".o", ".e")

                batch_error     = os.path.join(dir,      batch_error_base)
                batch_error_new = os.path.join(dir, id + batch_error_base)

                try:
                    os.rename(batch_error, batch_error_new)
                except FileNotFoundError:
                    pass

                break

        f.close()

    return

if __name__ == "__main__":
    if len(sys.argv) == 2:
        print("adding ids to batch files in dir: ", sys.argv[1])
        add_id_to_batch_outputs(dirname=sys.argv[1])
