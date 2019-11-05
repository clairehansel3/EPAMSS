#!/bin/sh
ln ../../scripts/run.py .
ln ../../scripts/simulation.py .
ln ../../epamss .
mkdir data
mkdir results
python3 run.py
echo "hash:"
md5sum --check data.md5
