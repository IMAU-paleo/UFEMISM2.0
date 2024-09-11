function analyse_component_tests_discretisation( ...
  foldername_component_tests_discretisation, foldername_automated_testing, do_print_figures)
% Analyse the results of all the discretisation component tests

analyse_component_tests_discretisation_map_deriv( ...
  [foldername_component_tests_discretisation '/mapping_derivatives'], foldername_automated_testing, do_print_figures)

end