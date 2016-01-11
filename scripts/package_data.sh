#!/bin/bash

tar -czvf data.tar.gz $(ls *.batch.* *.dat *.log 2>/dev/null)