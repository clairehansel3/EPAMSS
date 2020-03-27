#!/bin/bash
c++ -Wall -Wextra -std=c++17 -O3 -flto -ffast-math -DNDEBUG -I ../../source \
  ../../source/moliere.cxx ../../source/random.cxx moliere.cxx -o \
  moliere
./moliere
python3 moliere.py
