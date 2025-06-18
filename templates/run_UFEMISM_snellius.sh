#!/bin/bash
#Set job requirements
#SBATCH --time=05:00
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1

# Load modules for MPI, NetCDF, PETSc, and other required libraries
module load 2024
module load netCDF-Fortran/4.6.1-gompi-2024a
module load PETSc/3.22.0-foss-2024a
module load OpenMPI/5.0.3-GCC-13.3.0

# Execute the program
srun UFEMISM_program config-files/config_test
