#! /bin/csh -f

# Execute the program using -n N cores
mpiexec -n 2 UFEMISM_program config-files/config_test
