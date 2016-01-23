import platform

is_hyades = ("hyades" in platform.node()) or ("eudora" in platform.node())

if is_hyades:
    data_dir_default = "../cluster_parameter_study/"
else:
    data_dir_default = "../saved_runs/cluster_parameter_study/"