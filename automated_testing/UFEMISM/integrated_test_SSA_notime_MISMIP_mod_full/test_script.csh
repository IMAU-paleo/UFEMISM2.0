#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_SSA_notime_MISMIP_mod_full

rm -rf ${test_dir}/results*
mpiexec  -n 2  UFEMISM_program  ${test_dir}/config.cfg

python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}
