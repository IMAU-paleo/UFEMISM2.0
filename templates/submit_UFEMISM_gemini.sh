#!/bin/bash

# This is the script you run from the Gemini terminal
# to start a UFEMISM simulation.
# Doesn't need any arguments, those are listed below.
# The number of cores is set in run_UFEMISM_gemini.sh

# Load modules for MPI and PETSc
module load mpi/openmpi-x86_64
module load petsc/3.16.3
# NOTE: If you still get an 'mpiexec commmand not found'
# error, try loading the modules in your main terminal
# session (i.e. not within this submit script) before
# submiting your job. Good luck.

# Submit the run script with sbatch
qsub -cwd -m eas -V -q long.q ./run_UFEMISM_gemini.sh
