#!/bin/bash
#Set job requirements
#SBATCH --time=05:00
#SBATCH --partition=thin
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1

# Load modules for MPI, NetCDF, PETSc, and other required libraries
module load 2021
module load eb/4.5.2
module load foss/2021a
module load netCDF-Fortran/4.5.3-gompi-2021a
module load PETSc/3.15.1-foss-2021a
module load imkl/2021.2.0-iompi-2021a
module load OpenMPI/4.1.1-GCC-10.3.0

# Execute the program
srun UFEMISM_program config-files/config_test
