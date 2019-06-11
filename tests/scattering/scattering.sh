mkdir data
mkdir results
ln -s ../../scripts/simulation.py
make
python3 scattering.py
