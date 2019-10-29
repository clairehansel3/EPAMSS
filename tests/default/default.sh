#!/bin/sh
make
cp ../../scripts/run.py .
cp ../../scripts/simulation.py .
mkdir data
mkdir results
python3 run.py
echo "hash:"
md5sum --check data.md5
