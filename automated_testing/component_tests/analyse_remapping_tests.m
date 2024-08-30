function analyse_remapping_tests( foldername)

% foldername = '../../results_component_tests/remapping';
% addpath('../../tools/matlab/')

addpath('tools/matlab/')

do_print = false;

analyse_remapping_tests_grid_to_mesh( [foldername '/grid_to_mesh'], do_print);
analyse_remapping_tests_mesh_to_grid( [foldername '/mesh_to_grid'], do_print);

end