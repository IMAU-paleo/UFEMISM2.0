function analyse_component_tests( foldername_component_tests, foldername_automated_testing, do_print_figures)
% Analyse the results of all the component tests

analyse_component_tests_discretisation( ...
  [foldername_component_tests '/discretisation'], foldername_automated_testing, do_print_figures);

analyse_component_tests_remapping( ...
  [foldername_component_tests '/remapping'], foldername_automated_testing, do_print_figures);

end