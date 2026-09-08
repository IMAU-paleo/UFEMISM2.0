#! /bin/csh -f

set test_dir = automated_testing/UFEMISM/integrated_test_SSA_notime_MISMIP_mod_full

set solvers = (SSA SSA_FEM_PETSc)

rm -rf $test_dir/results*
mkdir $test_dir/results

foreach solver ($solvers)

  set exp_output_dir = $test_dir/results_${solver}

  rm -rf $exp_output_dir

  mpiexec -n 2 UFEMISM_program $test_dir/config_${solver}.cfg

  mv $exp_output_dir/main_output_ANT_00001.nc ${test_dir}/results/results_ISMIP_HOM_${solver}_mesh.nc
  mv $exp_output_dir/main_output_ANT_grid.nc  ${test_dir}/results/results_ISMIP_HOM_${solver}_grid.nc
  mv $exp_output_dir/checksum_logfile.txt     ${test_dir}/results/checksum_logfile_${solver}.txt

  rm -rf $exp_output_dir

end

python3 automated_testing/reduce_all_netcdfs_in_folder_to_checksum.py ${test_dir}
