#!/bin/sh
mkdir data
mkdir results
cp ../../scripts/simulation.py .
make clean
make
python3 scattering.py run
python3 scattering.py analyze
