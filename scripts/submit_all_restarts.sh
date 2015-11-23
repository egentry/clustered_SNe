#!/bin/bash

for f in "$( dirname "${BASH_SOURCE[0]}")/restart.batch.*"; do qsub $f; done
