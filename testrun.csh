#! /bin/csh -f

./compile_all_mac.csh

rm -rf results_UFEMISM_test
rm -rf results_UFEMISM_stresstest*
rm -rf results_2023*
rm -rf *.nc

mpiexec  -n 2  UFEMISM_program  config-files/config_test_01.cfg
#mpiexec  -n 2  UFEMISM_program  config-files/config_stresstest_01.cfg