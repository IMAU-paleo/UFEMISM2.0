#!/bin/sh

# Load necessary modules
module load 2021
module load eb/4.5.2
eblocalinstall PETSc-3.15.1-foss-2021a.eb
module load foss/2021a
module load netCDF-Fortran/4.5.3-gompi-2021a
module load PETSc/3.15.1-foss-2021a
module load imkl/2021.2.0-iompi-2021a
module load OpenMPI/4.1.1-GCC-10.3.0

# Go to src/, make, and come back
cd src
make all
cd ..

# Update program
rm -f UFEMISM_program
mv src/UFEMISM_program .
