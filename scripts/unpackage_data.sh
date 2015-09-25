#!/bin/bash


tar -xzvf data.tar.gz


# dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# echo $(dir)
python $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../analysis/parameter_study_file_structure.py $(pwd)