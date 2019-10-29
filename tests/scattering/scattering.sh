#!/bin/sh
mkdir data
mkdir results
cp ../../scripts/simulation.py .
cp ../../scripts/run_epamss.sh .
make
cp epamss ../../epamss
python3 scattering.py run
