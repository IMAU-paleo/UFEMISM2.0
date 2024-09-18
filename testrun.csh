#! /bin/csh -f

./compile_all_mac.csh

rm -rf automated_testing/integrated_tests/idealised/Halfar_dome/Halfar_5km/results

mpiexec  -n 2  UFEMISM_program  automated_testing/integrated_tests/idealised/Halfar_dome/Halfar_5km/config.cfg