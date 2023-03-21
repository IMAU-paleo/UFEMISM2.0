#! /bin/csh -f

./compile_all_mac.csh

rm -f testfile.nc
rm -f mesh_01.txt
rm -f mesh_02.txt


mpiexec  -n 2   UFEMISM_program