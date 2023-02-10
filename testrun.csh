#! /bin/csh -f

./compile_all_mac.csh

rm -f testfile.nc


mpiexec  -n 2   UFEMISM_program