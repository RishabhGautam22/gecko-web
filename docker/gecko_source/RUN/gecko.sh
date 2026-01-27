#!/bin/bash

# Exit on error
set -e

# Increase stack size for large mechanisms
ulimit -s unlimited || echo "Warning: Could not set ulimit -s unlimited"

cp ../INPUT/gecko.nml ./
cp ../INPUT/cheminput.dat ./
cp ../OBJ/cm ./

./cm

rm cm gecko.nml cheminput.dat

