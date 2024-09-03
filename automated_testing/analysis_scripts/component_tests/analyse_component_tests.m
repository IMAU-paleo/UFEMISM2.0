function analyse_component_tests( foldername_component_tests, do_print_figures)
% Analyse the results of all the component tests

% analyse_component_tests_discretisation( ...
%   [foldername_component_tests '/discretisation'], do_print_figures);

analyse_component_tests_remapping( ...
  [foldername_component_tests '/remapping'], do_print_figures);

end