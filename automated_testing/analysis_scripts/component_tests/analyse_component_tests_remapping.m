function analyse_component_tests_remapping( foldername, foldername_automated_testing, do_print_figures)
% Analyse the results of all the remapping component tests

analyse_component_tests_remapping_grid_to_mesh( ...
  [foldername '/grid_to_mesh'], foldername_automated_testing, do_print_figures);

analyse_component_tests_remapping_mesh_to_grid( ...
  [foldername '/mesh_to_grid'], foldername_automated_testing, do_print_figures);

analyse_component_tests_remapping_mesh_to_mesh( ...
  [foldername '/mesh_to_mesh'], foldername_automated_testing, do_print_figures);

end