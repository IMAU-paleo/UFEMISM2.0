#!/bin/sh

# Load necessary modules
module load mpi/openmpi-x86_64
module load petsc/3.16.3

# Go to src/, make, and come back
cd src
make clean
make all
cd ..

# Update program
rm -f UFEMISM_program
mv src/UFEMISM_program .
