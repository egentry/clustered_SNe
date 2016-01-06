#!/bin/bash

tar -xzvf data.tar.gz

python $( dirname "${BASH_SOURCE[0]}")/../analysis/add_id_to_batch_outputs.py $(pwd)
