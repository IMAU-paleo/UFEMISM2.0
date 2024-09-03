function analyse_component_tests_remapping( foldername, do_print_figures)
% Analyse the results of all the remapping component tests

analyse_component_tests_remapping_grid_to_mesh( ...
  [foldername '/grid_to_mesh'], do_print_figures);

analyse_component_tests_remapping_mesh_to_grid( ...
  [foldername '/mesh_to_grid'], do_print_figures);

end