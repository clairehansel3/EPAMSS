#!/bin/sh
mkdir data
mkdir results
#cp ../../scripts/simulation.py .
cp ../../scripts/run.py .
make
python3 scattering.py run
python3 scattering.py analyze
